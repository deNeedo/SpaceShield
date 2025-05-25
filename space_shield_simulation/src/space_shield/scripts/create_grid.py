import geopandas as gpd
from shapely.geometry import Polygon, box
import math
import os # For path joining
import folium # Added for visualization

def create_grid_for_polygon(polygon: Polygon, cell_size_degrees: float):
    """
    Creates a grid of square cells covering a given polygon.

    Args:
        polygon: The shapely Polygon to create a grid for.
        cell_size_degrees: The width and height of each grid cell in degrees.

    Returns:
        A list of shapely Polygons representing the grid cells that intersect
        the input polygon, clipped to the polygon's boundary.
    """
    minx, miny, maxx, maxy = polygon.bounds
    grid_cells = []

    # Iterate over the bounding box of the polygon with the given cell size
    x = minx
    while x < maxx:
        y = miny
        while y < maxy:
            # Create a square cell
            cell = box(x, y, x + cell_size_degrees, y + cell_size_degrees)
            
            # Check if the cell intersects with the input polygon
            if polygon.intersects(cell):
                # If it intersects, add the intersection part to the list
                # This clips the cell to the polygon's boundary
                intersection = polygon.intersection(cell)
                if not intersection.is_empty:
                    # Ensure we are adding actual polygons, not other geometry types
                    if intersection.geom_type == 'Polygon' or intersection.geom_type == 'MultiPolygon':
                        grid_cells.append(intersection)
            y += cell_size_degrees
        x += cell_size_degrees
    
    return grid_cells

def visualize_grid_on_map(original_aoi_gdf: gpd.GeoDataFrame, grid_gdf: gpd.GeoDataFrame, output_html_path: str):
    """
    Creates an HTML map visualizing the original AOI and the generated grid.
    Args:
        original_aoi_gdf: GeoDataFrame of the original Area of Interest.
        grid_gdf: GeoDataFrame of the generated grid cells.
        output_html_path: Path to save the HTML map.
    """
    print(f"Attempting to visualize grid on map: {output_html_path}")
    if original_aoi_gdf.empty:
        print("Original AOI GeoDataFrame is empty, cannot determine map center.")
        # Fallback map center if AOI is empty
        map_center = [0,0]
        zoom_start = 2
    else:
        # Calculate a suitable map center and zoom
        # Use the centroid of the union of all original geometries
        map_center_geom = original_aoi_gdf.unary_union.centroid
        map_center = [map_center_geom.y, map_center_geom.x]
        zoom_start = 12 # Adjust as needed

    m = folium.Map(location=map_center, zoom_start=zoom_start, tiles="OpenStreetMap")

    # Add original AOI polygons to the map
    if not original_aoi_gdf.empty:
        folium.GeoJson(
            original_aoi_gdf,
            name='Original AOI',
            style_function=lambda x: {'fillColor': 'blue', 'color': 'blue', 'weight': 2, 'fillOpacity': 0.1}
        ).add_to(m)

    # Add grid cells to the map
    if not grid_gdf.empty:
        folium.GeoJson(
            grid_gdf,
            name='Generated Grid',
            style_function=lambda x: {'fillColor': 'red', 'color': 'black', 'weight': 1, 'fillOpacity': 0.4},
            tooltip=folium.features.GeoJsonTooltip(fields=['cell_id'], aliases=['Cell ID:'], localize=True)
        ).add_to(m)
    else:
        print("Grid GeoDataFrame is empty, nothing to add to map for grid.")

    folium.LayerControl().add_to(m)

    try:
        m.save(output_html_path)
        print(f"Successfully saved map visualization to {output_html_path}")
    except Exception as e:
        print(f"Error saving map HTML: {e}")

def main():
    # --- Configuration ---
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    
    input_geojson_path = os.path.join(project_root, "data", "simulation_inputs", "input_flood.geojson")
    output_grid_geojson_path = os.path.join(project_root, "data", "simulation_outputs", "flood_area_grid.geojson")
    output_map_html_path = os.path.join(project_root, "data", "simulation_outputs", "flood_area_grid_map.html") # Path for the map
    
    cell_size_degrees = 0.02

    # --- Load Input GeoJSON ---
    try:
        aoi_gdf = gpd.read_file(input_geojson_path)
        print(f"Successfully loaded {input_geojson_path}")
        print(f"CRS of input GeoJSON: {aoi_gdf.crs}")
    except Exception as e:
        print(f"Error loading GeoJSON file: {e}")
        return

    if aoi_gdf.empty:
        print("Input GeoJSON is empty.")
        return

    # --- Generate Grid Cells ---
    all_grid_cells = []
    for index, row in aoi_gdf.iterrows():
        polygon_geometry = row['geometry']
        if isinstance(polygon_geometry, Polygon):
            # print(f"Processing polygon {index+1}...")
            grid_cells_for_this_polygon = create_grid_for_polygon(polygon_geometry, cell_size_degrees)
            all_grid_cells.extend(grid_cells_for_this_polygon)
        else:
            if hasattr(polygon_geometry, 'geoms'): 
                # print(f"Processing MultiPolygon {index+1}...")
                for poly_part in polygon_geometry.geoms:
                    if isinstance(poly_part, Polygon):
                        grid_cells_for_this_polygon = create_grid_for_polygon(poly_part, cell_size_degrees)
                        all_grid_cells.extend(grid_cells_for_this_polygon)
                    # else:
                        # print(f"Skipping non-Polygon part within MultiPolygon at index {index}: {type(poly_part)}")
            # else:
                # print(f"Skipping non-Polygon geometry at index {index}: {type(polygon_geometry)}")

    if not all_grid_cells:
        print("No grid cells generated. Check if polygons exist and cell size is appropriate.")
        # Still try to visualize the original AOI if it exists
        if not aoi_gdf.empty:
            visualize_grid_on_map(aoi_gdf, gpd.GeoDataFrame(geometry=[], crs=aoi_gdf.crs), output_map_html_path)
        return

    # --- Create GeoDataFrame from Grid Cells ---
    grid_gdf = gpd.GeoDataFrame(geometry=all_grid_cells, crs=aoi_gdf.crs)
    grid_gdf['cell_id'] = range(len(grid_gdf))
    print(f"Generated {len(grid_gdf)} grid cells.")

    # --- Save Output Grid GeoJSON ---
    output_dir = os.path.dirname(output_grid_geojson_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    try:
        grid_gdf.to_file(output_grid_geojson_path, driver="GeoJSON")
        print(f"Successfully saved grid to {output_grid_geojson_path}")
    except Exception as e:
        print(f"Error saving grid GeoJSON: {e}")

    # --- Visualize on Map ---
    visualize_grid_on_map(aoi_gdf, grid_gdf, output_map_html_path)

if __name__ == "__main__":
    main() 