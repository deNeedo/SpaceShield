import os
import http.server
import socketserver
import webbrowser
import argparse
import json # Added for loading GeoJSON

import folium
from folium.plugins import Draw

# Define the directory for our local web server content and map filename
# Assumes this script is in space_shield_simulation/src/space_shield/
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")) # Points to space_shield_simulation
LOCAL_WEB_DIR = os.path.join(BASE_DIR, "local_web")
MAP_FILENAME = "map.html"
SERVER_PORT = 8000

DEFAULT_LAT = 50.372003383912485  # Ljubljana Latitude
DEFAULT_LON = 16.67052239355332  # Ljubljana Longitude

def generate_map_file(lat, lon, output_path, geojson_aoi_path=None):
    """Generates the folium map and saves it to output_path.
    Optionally loads and displays a GeoJSON AOI.
    """
    print(f"Generating map centered at Latitude: {lat}, Longitude: {lon}")
    if geojson_aoi_path:
        print(f"Attempting to load AOI from: {geojson_aoi_path}")
        
    folium_map = folium.Map(location=[lat, lon], zoom_start=10, tiles="OpenStreetMap")
    
    # Add various tile layers for user selection
    folium.TileLayer(
        tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr='Esri', name='Esri Satellite', overlay=False, control=True, show=False
    ).add_to(folium_map)
    folium.TileLayer(
        tiles='https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
        attr='Google', name='Google Satellite', overlay=False, control=True, show=True
    ).add_to(folium_map)
    folium.TileLayer(
        tiles='https://mt1.google.com/vt/lyrs=m&x={x}&y={y}&z={z}',
        attr='Google', name='Google Maps', overlay=False, control=True, show=False
    ).add_to(folium_map)
    folium.TileLayer(
        tiles='https://mt1.google.com/vt/lyrs=p&x={x}&y={y}&z={z}',
        attr='Google', name='Google Terrain', overlay=False, control=True, show=False
    ).add_to(folium_map)

    # Add drawing tools
    draw_plugin = Draw(
        export=True, filename='drawn_aoi.geojson', position='topleft',
        draw_options={
            'polyline': False, 'rectangle': True, 'circle': False,
            'marker': False, 'circlemarker': False,
            'polygon': {'allowIntersection': False}
        },
        edit_options={'edit': True, 'remove': True}
    ).add_to(folium_map)
    
    # Load and display AOI if path is provided and file exists
    if geojson_aoi_path and os.path.exists(geojson_aoi_path):
        try:
            with open(geojson_aoi_path, 'r') as f:
                aoi_data = json.load(f)
            folium.GeoJson(
                aoi_data,
                name='Loaded AOI',
                style_function=lambda x: {'fillColor': 'orange', 'color': 'orange', 'weight': 2.5, 'fillOpacity': 0.4}
            ).add_to(folium_map)
            print(f"Successfully loaded and added AOI from {geojson_aoi_path} to map.")
        except Exception as e_aoi:
            print(f"[WARNING] Could not load or display AOI from {geojson_aoi_path}: {e_aoi}")
    elif geojson_aoi_path:
        print(f"[WARNING] AOI file not found: {geojson_aoi_path}")

    folium.LayerControl().add_to(folium_map)

    try:
        # Ensure the local_web directory exists
        if not os.path.exists(LOCAL_WEB_DIR):
            os.makedirs(LOCAL_WEB_DIR)
            print(f"Created directory: {LOCAL_WEB_DIR}")
            
        folium_map.save(output_path)
        print(f"Map file saved to: {output_path}")
    except Exception as e:
        print(f"[ERROR] Could not save map file: {e}")
        import traceback
        traceback.print_exc()

def start_server(directory, port):
    """Starts a simple HTTP server in the specified directory and port."""
    class CustomHandler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=directory, **kwargs)

    try:
        with socketserver.TCPServer(("localhost", port), CustomHandler) as httpd:
            print(f"HTTP server started on http://localhost:{port}/")
            print(f"Serving files from: {directory}")
            map_url = f"http://localhost:{port}/{MAP_FILENAME}"
            print(f"Map should be available at: {map_url}")
            
            try:
                webbrowser.open_new_tab(map_url)
                print(f"Attempted to open map in web browser.")
            except Exception as e_wb:
                print(f"Could not automatically open web browser: {e_wb}. Please open manually.")

            print("Press Ctrl+C to stop the server.")
            httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopping...")
    except Exception as e:
        print(f"[ERROR] Could not start server: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a Folium map and serve it locally.")
    parser.add_argument("--lat", type=float, default=DEFAULT_LAT, help=f"Center latitude for the map (default: {DEFAULT_LAT})")
    parser.add_argument("--lon", type=float, default=DEFAULT_LON, help=f"Center longitude for the map (default: {DEFAULT_LON})")
    parser.add_argument("--port", type=int, default=SERVER_PORT, help=f"Port for the local HTTP server (default: {SERVER_PORT})")
    parser.add_argument("--aoi", type=str, help="Path to a GeoJSON file to display as an Area of Interest.")
    
    args = parser.parse_args()

    map_file_full_path = os.path.join(LOCAL_WEB_DIR, MAP_FILENAME)
    
    # Generate the map first
    generate_map_file(args.lat, args.lon, map_file_full_path, geojson_aoi_path=args.aoi)
    
    # Then start the server
    if os.path.exists(map_file_full_path):
        start_server(LOCAL_WEB_DIR, args.port)
    else:
        print(f"[ERROR] Map file {map_file_full_path} was not generated. Cannot start server.") 