# Geospatial utilities for handling coordinates, polygons, etc.

import json
from shapely.geometry import Polygon, Point

def load_aoi_from_geojson(filepath: str) -> Polygon:
    """Loads an Area of Interest (AOI) polygon from a GeoJSON file."""
    try:
        with open(filepath, 'r') as f:
            geojson_data = json.load(f)
        
        # Assuming the first feature is the AOI polygon
        # This might need to be more robust depending on GeoJSON structure
        coordinates = geojson_data['features'][0]['geometry']['coordinates'][0]
        return Polygon(coordinates)
    except Exception as e:
        print(f"Error loading AOI from {filepath}: {e}")
        return None

def create_rectangular_aoi(min_x: float, min_y: float, max_x: float, max_y: float) -> Polygon:
    """Creates a rectangular AOI polygon."""
    return Polygon([
        (min_x, min_y),
        (min_x, max_y),
        (max_x, max_y),
        (max_x, min_y),
        (min_x, min_y)
    ])

# Add other utility functions as needed:
# - Coordinate transformations (if dealing with lat/lon vs. local cartesian)
# - Distance calculations (e.g., haversine if using lat/lon)
# - Point in polygon checks (Shapely does this with polygon.contains(point)) 