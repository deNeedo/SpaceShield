import geopandas as gpd
from shapely.geometry import Point, Polygon, box, LineString
from shapely.ops import unary_union
import numpy as np
import folium
from folium.plugins import TimestampedGeoJson # Added for animation
import math
import os
import time
from datetime import datetime, timedelta # Added for animation timing

from python_tsp.heuristics import solve_tsp_simulated_annealing
from python_tsp.distances import great_circle_distance_matrix
import pulp

# --- Configuration & Parameters ---
DRONE_COVERAGE_RADIUS_KM = 10.0 # Drone coverage radius in kilometers
DRONE_COVERAGE_RADIUS_KM_actual = 15.0

# --- Global Helper for Geopy (for accurate distance calculation) ---
_geopy_available = False
try:
    from geopy.distance import great_circle
    _geopy_available = True
    print("INFO: geopy.distance found. Using for accurate path distance calculations in animation.")
except ImportError:
    print("WARNING: geopy.distance module not found. Animation will use less accurate fallback for path distance calculations. Drone spacing and speed may not be precise. Consider `pip install geopy`.")

# --- Helper function to create a circular polygon from a center and radius ---
def create_circle_geo(center_point: Point, radius_km: float, segments=32) -> Polygon:
    coords = []
    center_lon, center_lat = center_point.x, center_point.y
    radius_deg_lat = radius_km / 111.1
    
    effective_lat_for_cos = center_lat
    if center_lat >= 89.99: effective_lat_for_cos = 89.0
    elif center_lat <= -89.99: effective_lat_for_cos = -89.0
    
    cos_lat = math.cos(math.radians(effective_lat_for_cos))
    if cos_lat < 0.000001:
        radius_deg_lon = radius_deg_lat * 10 
    else:
        radius_deg_lon = radius_km / (111.1 * cos_lat)

    for i in range(segments + 1):
        angle_rad = (i / segments) * 2 * math.pi
        dx_deg = radius_deg_lon * math.cos(angle_rad) 
        dy_deg = radius_deg_lat * math.sin(angle_rad) 
        coords.append((center_lon + dx_deg, center_lat + dy_deg))
    return Polygon(coords)

# --- Helper Functions for Animation Path Calculations ---
def calculate_path_segments_and_total_length(path_coords: list) -> tuple[list, float]:
    """Calculates lengths of each segment and total length of the path in km."""
    segment_lengths_km = []
    total_length_km = 0.0
    if len(path_coords) < 2:
        return [], 0.0

    for i in range(len(path_coords) - 1):
        p1_lon, p1_lat = path_coords[i]
        p2_lon, p2_lat = path_coords[i+1]
        
        segment_len = 0.0
        if _geopy_available:
            segment_len = great_circle((p1_lat, p1_lon), (p2_lat, p2_lon)).km
        else: 
            dy_km = (p2_lat - p1_lat) * 111.1
            dx_km = (p2_lon - p1_lon) * 111.1 * math.cos(math.radians((p1_lat + p2_lat) / 2.0))
            segment_len = math.sqrt(dx_km**2 + dy_km**2)
        
        segment_lengths_km.append(segment_len)
        total_length_km += segment_len
    return segment_lengths_km, total_length_km

def get_point_at_distance_on_path(
    path_coords: list, 
    distance_km: float, 
    segment_lengths_km: list, 
    total_path_length_km: float
) -> Point | None:
    """
    Finds a point along a path (list of (lon, lat) tuples) at a specific distance in km.
    Uses pre-calculated segment lengths and total length. Handles path wrapping.
    """
    if not path_coords:
        print("Error: Path coordinates are empty for animation point calculation.")
        return None
    if total_path_length_km == 0.0: # Single point path or all points identical
        return Point(path_coords[0])

    target_dist_km = distance_km % total_path_length_km
    if target_dist_km < 0: 
        target_dist_km += total_path_length_km
    
    epsilon = 1e-9 
    if target_dist_km < epsilon or abs(target_dist_km - total_path_length_km) < epsilon : # At start or very end
         return Point(path_coords[0]) # Path is closed, so end is start

    current_distance_sum_km = 0.0
    for i in range(len(segment_lengths_km)):
        segment_len = segment_lengths_km[i]
        if current_distance_sum_km + segment_len >= target_dist_km - epsilon:
            if segment_len < epsilon: 
                return Point(path_coords[i]) 

            remaining_dist_on_segment = target_dist_km - current_distance_sum_km
            remaining_dist_on_segment = max(0.0, min(remaining_dist_on_segment, segment_len))
            
            fraction = remaining_dist_on_segment / segment_len if segment_len > epsilon else 0.0
            
            start_node_lon, start_node_lat = path_coords[i]
            end_node_lon, end_node_lat = path_coords[i+1]

            interp_lon = start_node_lon + fraction * (end_node_lon - start_node_lon)
            interp_lat = start_node_lat + fraction * (end_node_lat - start_node_lat)
            return Point(interp_lon, interp_lat)
        current_distance_sum_km += segment_len
    
    return Point(path_coords[-1]) # Fallback for end of path (same as start for closed path)

def main():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) if '__file__' in globals() else os.getcwd()
    data_dir = os.path.join(project_root, "data", "simulation_outputs")
    os.makedirs(data_dir, exist_ok=True) # Ensure output directory exists

    grid_geojson_path = os.path.join(data_dir, "flood_area_grid.geojson")
    output_circles_geojson_path = os.path.join(data_dir, "coverage_circles.geojson")
    output_path_geojson_path = os.path.join(data_dir, "coverage_path.geojson")
    output_map_html_path = os.path.join(data_dir, "coverage_circles_map.html")
    animated_map_path_html = os.path.join(data_dir, "animated_drones_map.html") # Path for animated map

    # --- 1. Load Grid Cells ---
    try:
        # Create a dummy grid_gdf if file doesn't exist for testing purposes
        if not os.path.exists(grid_geojson_path):
            print(f"Warning: Grid file {grid_geojson_path} not found. Creating a dummy grid for demonstration.")
            # Example: A 2x2 grid around a sample point
            base_lon, base_lat = 16.0, 50.0
            points = [
                Point(base_lon, base_lat), Point(base_lon + 0.1, base_lat),
                Point(base_lon, base_lat + 0.1), Point(base_lon + 0.1, base_lat + 0.1),
                Point(base_lon + 0.2, base_lat + 0.2) # Add more to ensure enough cells for TSP
            ]
            grid_gdf = gpd.GeoDataFrame(
                {'cell_id': range(len(points))},
                geometry=[p.buffer(0.05, cap_style=3) for p in points], # Create square cells
                crs="EPSG:4326"
            )
            grid_gdf.to_file(grid_geojson_path, driver="GeoJSON")
            print(f"Dummy grid created and saved to {grid_geojson_path}.")
        
        grid_gdf = gpd.read_file(grid_geojson_path)
        if grid_gdf.empty:
            print(f"Grid file {grid_geojson_path} is empty. Exiting.")
            return
        if 'cell_id' not in grid_gdf.columns:
            grid_gdf['cell_id'] = range(len(grid_gdf))
        print(f"Loaded {len(grid_gdf)} grid cells from {grid_geojson_path}")
    except Exception as e:
        print(f"Error loading or creating grid GeoJSON: {e}")
        return

    all_cell_ids = grid_gdf['cell_id'].tolist()
    
    # --- 2. Generate All Potential Candidate Circles and Determine Coverage ---
    print("Generating all potential candidate circles and their coverage...")
    candidate_circles_info = []
    for idx, cell_row in grid_gdf.iterrows():
        potential_center = cell_row.geometry.centroid
        potential_circle_poly = create_circle_geo(potential_center, DRONE_COVERAGE_RADIUS_KM)
        covered_cell_ids_by_this_candidate = set()
        for _, check_cell_row in grid_gdf.iterrows(): # Use unique var name for inner loop
            if potential_circle_poly.covers(check_cell_row.geometry.centroid): # Check if centroid is covered
                covered_cell_ids_by_this_candidate.add(check_cell_row['cell_id'])
        
        if covered_cell_ids_by_this_candidate:
            candidate_circles_info.append({
                "center": potential_center, "polygon": potential_circle_poly,
                "covers": covered_cell_ids_by_this_candidate, "id": f"candidate_circle_{idx}"
            })
    
    if not candidate_circles_info:
        print("No potential candidate circles could be generated or none cover any cells. Exiting.")
        return
    print(f"Generated {len(candidate_circles_info)} potential candidate circles.")

    # --- 3. Solve Set Cover Problem using Integer Linear Programming (PuLP) ---
    print("Setting up and solving the Set Cover problem using PuLP...")
    start_time_pulp = time.time()
    prob = pulp.LpProblem("DroneSetCover", pulp.LpMinimize)
    circle_vars = pulp.LpVariable.dicts("Circle", [info["id"] for info in candidate_circles_info], cat=pulp.LpBinary)
    prob += pulp.lpSum([circle_vars[info["id"]] for info in candidate_circles_info]), "TotalCircles"

    for cell_id in all_cell_ids:
        covering_circles_for_this_cell = [circle_vars[info["id"]] for info in candidate_circles_info if cell_id in info["covers"]]
        if covering_circles_for_this_cell:
             prob += pulp.lpSum(covering_circles_for_this_cell) >= 1, f"Cell_{cell_id}_Covered"
        else:
             print(f"Warning: Cell ID {cell_id} cannot be covered by any potential circle. Check grid or radius.")
    
    prob.solve(pulp.PULP_CBC_CMD(msg=0)) # Suppress solver messages
    pulp_duration = time.time() - start_time_pulp
    print(f"PuLP solver status: {pulp.LpStatus[prob.status]}. Solved in {pulp_duration:.2f} seconds.")

    if prob.status != pulp.LpStatusOptimal:
        print("Optimal solution not found by PuLP. Exiting.")
        # (Error details as before)
        return

    chosen_circles_centers = []
    chosen_circles_polygons = []
    for info in candidate_circles_info:
        if pulp.value(circle_vars[info["id"]]) == 1:
            chosen_circles_centers.append(info["center"])
            chosen_circles_polygons.append(info["polygon"])
    print(f"\nILP algorithm finished. Optimal number of circles found: {len(chosen_circles_polygons)}.")
    
    # --- 4. Save and Visualize (Static Map) --- 
    if not chosen_circles_polygons:
        print("No coverage circles were selected by the ILP solver.")
        return

    circles_gdf = gpd.GeoDataFrame(geometry=chosen_circles_polygons, crs=grid_gdf.crs)
    circles_gdf['circle_id'] = range(len(circles_gdf))
    try:
        circles_gdf.to_file(output_circles_geojson_path, driver="GeoJSON")
        print(f"Successfully saved coverage circles to {output_circles_geojson_path}")
    except Exception as e:
        print(f"Error saving coverage circles GeoJSON: {e}")

    path_gdf = None
    path_line = None # Define path_line here to be accessible for animation
    if len(chosen_circles_centers) > 1:
        coords_array = np.array([[pt.y, pt.x] for pt in chosen_circles_centers])
        distance_matrix_val = great_circle_distance_matrix(coords_array)
        permutation, _ = solve_tsp_simulated_annealing(distance_matrix_val)
        optimized_centers = [chosen_circles_centers[i] for i in permutation]
        
        line_coords = [(pt.x, pt.y) for pt in optimized_centers]
        if len(line_coords) > 1: 
            line_coords.append(line_coords[0]) 
            path_line = LineString(line_coords) # Assign to path_line
            path_gdf = gpd.GeoDataFrame(geometry=[path_line], crs=grid_gdf.crs)
            path_gdf['path_id'] = ['coverage_path_01']
            try:
                path_gdf.to_file(output_path_geojson_path, driver="GeoJSON")
                print(f"Successfully saved coverage path to {output_path_geojson_path}")
            except Exception as e:
                print(f"Error saving coverage path GeoJSON: {e}")
    elif len(chosen_circles_centers) == 1:
        print("Only one circle center found, no path to create.")
    else:
        print("No circle centers to create a path from.")

    # Create Static Folium map
    map_center_geom = grid_gdf.unary_union.centroid if not grid_gdf.empty else Point(0,0).buffer(1).centroid
    map_center = [map_center_geom.y, map_center_geom.x]
    m_static = folium.Map(location=map_center, zoom_start=10, tiles="OpenStreetMap")
    folium.GeoJson(grid_gdf, name='AOI Grid Cells', style_function=lambda x: {'fillColor': 'gray', 'color': 'black', 'weight':0.5, 'fillOpacity':0.1}).add_to(m_static)
    folium.GeoJson(circles_gdf, name='Coverage Circles', style_function=lambda x: {'fillColor': 'red', 'color': 'red', 'weight':1, 'fillOpacity':0.3}, tooltip=folium.features.GeoJsonTooltip(fields=['circle_id'], aliases=['Circle ID:'])).add_to(m_static)
    if path_gdf is not None and not path_gdf.empty:
        folium.GeoJson(path_gdf, name='Coverage Path', style_function=lambda x: {'color': 'blue', 'weight': 3, 'opacity': 0.7}).add_to(m_static)
    folium.LayerControl().add_to(m_static)
    try:
        m_static.save(output_map_html_path)
        print(f"Successfully saved static coverage map to {output_map_html_path}")
    except Exception as e:
        print(f"Error saving static map HTML: {e}")

    # --- 5. Add Drone Animation ---
    if path_line and len(chosen_circles_centers) > 1 : # Check if path_line was created
        print("\nPreparing drone animation...")
        
        NUM_DRONES = 3
        DRONE_SPACING_KM = 15.0
        DRONE_ANIMATION_SPEED_KMH = 120.0 
        SIMULATION_DURATION_MINUTES = 30 
        FRAMES_PER_SIM_SECOND = 1 

        path_coords_list = list(path_line.coords)
        path_segment_lengths_km, current_total_path_length_km = calculate_path_segments_and_total_length(path_coords_list)

        if current_total_path_length_km == 0.0 and len(path_coords_list) > 1:
            print("Warning: Calculated path length for animation is 0. Drones may not move.")
            current_total_path_length_km = 1.0 # Avoid division by zero

        drone_speed_km_per_sec = DRONE_ANIMATION_SPEED_KMH / 3600.0
        total_simulation_seconds = int(SIMULATION_DURATION_MINUTES * 60)
        num_animation_frames = int(total_simulation_seconds * FRAMES_PER_SIM_SECOND)
        
        animation_features = []
        start_sim_datetime = datetime.utcnow()
        drone_colors = ['#FF0000', '#0000FF', '#008000'] # Red, Blue, Green

        print(f"Total path length for animation: {current_total_path_length_km:.2f} km")
        print(f"Animating {NUM_DRONES} drones for {SIMULATION_DURATION_MINUTES} minutes ({num_animation_frames} frames).")

        for frame_idx in range(num_animation_frames):
            sim_time_elapsed_s = frame_idx / FRAMES_PER_SIM_SECOND
            current_timestamp_dt = start_sim_datetime + timedelta(seconds=sim_time_elapsed_s)
            time_str_iso = current_timestamp_dt.isoformat() + "Z"

            distance_lead_drone_km = drone_speed_km_per_sec * sim_time_elapsed_s
            
            for drone_idx in range(NUM_DRONES):
                drone_effective_distance_km = distance_lead_drone_km - (drone_idx * DRONE_SPACING_KM)
                drone_center_point = get_point_at_distance_on_path(
                    path_coords_list, drone_effective_distance_km, 
                    path_segment_lengths_km, current_total_path_length_km
                )

                if drone_center_point is None: continue
                drone_color = drone_colors[drone_idx % len(drone_colors)]

                # Drone marker feature
                animation_features.append({
                    'type': 'Feature',
                    'geometry': drone_center_point.__geo_interface__,
                    'properties': {
                        'times': [time_str_iso],
                        'icon': 'circle',
                        'iconstyle': {
                            'fillColor': drone_color, 'color': '#000000', 'fillOpacity': 0.9, 
                            'weight': 1, 'radius': 6
                        },
                        'popup': f"Drone {drone_idx+1}<br>Time: {current_timestamp_dt.strftime('%H:%M:%S')}"
                    }
                })

                # Drone coverage circle feature
                coverage_poly = create_circle_geo(drone_center_point, DRONE_COVERAGE_RADIUS_KM_actual)
                animation_features.append({
                    'type': 'Feature',
                    'geometry': coverage_poly.__geo_interface__,
                    'properties': {
                        'times': [time_str_iso],
                        'style': {
                            'color': drone_color, 'weight': 1, 'opacity': 0.4,
                            'fillColor': drone_color, 'fillOpacity': 0.10,
                        },
                        'popup': f"Drone {drone_idx+1} Coverage"
                    }
                })
            
            if frame_idx % (60 * FRAMES_PER_SIM_SECOND) == 0:
                print(f"  Generated animation frame {frame_idx+1}/{num_animation_frames} (Sim time: {sim_time_elapsed_s/60.0:.1f} min)")

        animation_geojson_data = {'type': 'FeatureCollection', 'features': animation_features}
        m_animated = folium.Map(location=map_center, zoom_start=10, tiles="CartoDB positron")

        # Add static elements for context in animated map
        folium.GeoJson(grid_gdf, name='AOI Grid Cells', style_function=lambda x: {'fillColor': '# sosokosos', 'color': '#222222', 'weight':0.5, 'fillOpacity':0.05}).add_to(m_animated)
        if path_gdf is not None:
             folium.GeoJson(path_gdf, name='Coverage Path (Static)', style_function=lambda x: {'color': 'black', 'weight': 1.5, 'opacity': 0.6, 'dashArray': '5, 5'}).add_to(m_animated)
        
        TimestampedGeoJson(
            animation_geojson_data,
            period=f'PT{1.0/FRAMES_PER_SIM_SECOND:.0f}S', # Period between frames
            add_last_point=False, # Drones move, don't keep old positions
            auto_play=True, loop=True, loop_button=True,
            date_options='YYYY/MM/DD HH:mm:ss UTC',
            time_slider_drag_update=True,
            duration=f'PT{1}S',
            transition_time=int((1.0/FRAMES_PER_SIM_SECOND)*1000), # ms, display time for each frame's data
        ).add_to(m_animated)

        folium.LayerControl().add_to(m_animated)
        try:
            m_animated.save(animated_map_path_html)
            print(f"Successfully saved animated drone map to {animated_map_path_html}")
        except Exception as e:
            print(f"Error saving animated map HTML: {e}")
    else:
        print("\nSkipping drone animation: No valid TSP path found or too few nodes.")

if __name__ == "__main__":
    main()