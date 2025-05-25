# Data models for SpaceShield (e.g., Drone, AOI, BaseStation)

from dataclasses import dataclass, field
from typing import Tuple, List
from shapely.geometry import Point, Polygon
import folium
from folium.plugins import TimestampedGeoJson
import math
from datetime import datetime, timedelta

@dataclass
class BaseStation:
    id: str
    position: Point # Using Shapely Point for location
    # Add other relevant attributes like available communication channels etc.

@dataclass
class Drone:
    id: str
    # Configurable parameters
    speed_mps: float = 33.0
    coverage_radius_m: float = 15000.0
    effective_range_m: float = 15000.0
    max_flight_hours_h: float = 8.0

    # State variables
    current_position: Point = field(default_factory=lambda: Point(0,0))
    # current_battery_wh: float = field(init=False) # Replaced by flight_time_remaining_s
    flight_time_remaining_s: float = field(init=False)
    is_active: bool = True
    current_path: List[Tuple[float, float]] = field(default_factory=list)
    path_index: int = 0

    def __post_init__(self):
        self.flight_time_remaining_s = self.max_flight_hours_h * 3600.0

    def update_flight_time(self, time_elapsed_s: float):
        """Updates the remaining flight time and deactivates drone if time is up."""
        if not self.is_active:
            return
        self.flight_time_remaining_s -= time_elapsed_s
        if self.flight_time_remaining_s <= 0:
            self.is_active = False
            self.flight_time_remaining_s = 0
            print(f"Drone {self.id} max flight time reached.")

    def move_to(self, new_position: Point, time_elapsed_s: float):
        """Moves the drone and updates its flight time."""
        self.current_position = new_position
        self.update_flight_time(time_elapsed_s)

    def __repr__(self):
        return (
            f"Drone(id={self.id}, pos=({self.current_position.x:.2f}, {self.current_position.y:.2f}), "
            f"flight_time_rem={self.flight_time_remaining_s/3600.0:.2f}h, active={self.is_active})"
        )

@dataclass
class AreaOfInterest:
    id: str
    geometry: Polygon # Shapely Polygon
    # Add other attributes like priority, required coverage percentage, etc.

# You might also want models for:
# - SimulationState
# - PathSegment
# - CommunicationLink 

# --- Helper function to create coordinates for a circle polygon ---
def create_circle_polygon_coords(center_lat, center_lon, radius_km, segments=32):
    """
    Generates approximate geographic coordinates for a circle polygon.
    Note: This is an approximation. For high accuracy, especially for large radii
    or near poles, use a geospatial library.
    """
    coords = []
    # Approximate conversion: 1 degree latitude ~ 111.1 km
    # 1 degree longitude ~ 111.1 km * cos(latitude)
    radius_deg_lat = radius_km / 111.1
    radius_deg_lon = radius_km / (111.1 * math.cos(math.radians(center_lat)))

    for i in range(segments + 1):
        angle = (i / segments) * 2 * math.pi
        dy = radius_deg_lat * math.sin(angle)
        dx = radius_deg_lon * math.cos(angle)
        coords.append([center_lon + dx, center_lat + dy])
    return [coords] # GeoJSON Polygon coordinates are a list of lists

# --- Simulation Parameters ---
map_center_lat = 50.372003383912485
map_center_lon = 16.67052239355332
map_zoom = 9

plane_path_center_lat = 50.45
plane_path_center_lon = 16.8
plane_path_radius_deg = 0.20  # Radius of the plane's circular flight path in degrees
coverage_radius_km = 10.0

num_timesteps = 60*60  # Set for one update per second over 60 minutes
total_duration_minutes = 60  # 1 hour

# --- Generate GeoJSON features for TimestampedGeoJson ---
features = []
start_time = datetime.utcnow()

for i in range(num_timesteps):
    # Calculate plane's position on its circular path
    angle = (i / num_timesteps) * 2 * math.pi
    plane_lat = plane_path_center_lat + plane_path_radius_deg * math.sin(angle)
    plane_lon = plane_path_center_lon + plane_path_radius_deg * math.cos(angle)

    # Calculate timestamp for this specific point in time
    current_time_for_step = start_time + timedelta(seconds=(i * (total_duration_minutes * 60) / num_timesteps))
    time_str_for_step = current_time_for_step.isoformat() + "Z" # ISO format with Z for UTC

    # Create coverage polygon
    coverage_coords = create_circle_polygon_coords(plane_lat, plane_lon, coverage_radius_km)

    features.append({
        'type': 'Feature',
        'geometry': {   
            'type': 'Polygon',
            'coordinates': coverage_coords,
        },
        'properties': {
            'times': [time_str_for_step],
            'style': { # Style for a single, non-persistent circle
                'color': 'blue',
                'weight': 2,
                'opacity': 0.7,
                'fillColor': 'blue',
                'fillOpacity': 0.2,
            },
            'popup': f"Coverage at {current_time_for_step.strftime('%Y-%m-%d %H:%M:%S UTC')}<br>Plane Lat: {plane_lat:.4f}, Lon: {plane_lon:.4f}"
        }
    })

geojson_data = {
    'type': 'FeatureCollection',
    'features': features,
}

# --- Create Folium Map ---
# Use the initial map center and zoom from your existing map.html if desired
m = folium.Map(location=[map_center_lat, map_center_lon], zoom_start=map_zoom)

# Add some tile layers like in your original map
folium.TileLayer('openstreetmap', name='OpenStreetMap').add_to(m)
folium.TileLayer(
    tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
    attr='Esri',
    name='Esri Satellite',
    overlay=False,
    control=True
).add_to(m)
folium.TileLayer(
    tiles='https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
    attr='Google',
    name='Google Satellite',
    overlay=False,
    control=True
).add_to(m)


# --- Add TimestampedGeoJson to the map ---
TimestampedGeoJson(
    geojson_data,
    period='PT{:.0f}S'.format((total_duration_minutes * 60) / num_timesteps), # Period between frames
    add_last_point=False,
    auto_play=True,
    loop=True,
    loop_button=True,
    date_options='YYYY/MM/DD HH:mm:ss UTC',
    time_slider_drag_update=True,
    duration='PT{:.0f}M'.format((total_duration_minutes * 60) / num_timesteps), # Total duration of the animation
    # Transition time appropriate for non-persistent features
    transition_time = int(((total_duration_minutes * 60) / num_timesteps)) if num_timesteps > 0 else 200
).add_to(m)

# Add Layer Control to switch between base maps
folium.LayerControl().add_to(m)

# --- Save the map to an HTML file ---
output_filename = "animated_map.html"
m.save(output_filename)

print(f"Animated map saved to {output_filename}")
print(f"Plane path center: Lat={plane_path_center_lat}, Lon={plane_path_center_lon}")
print(f"Coverage radius: {coverage_radius_km} km")
print(f"Number of animation steps: {num_timesteps}") 