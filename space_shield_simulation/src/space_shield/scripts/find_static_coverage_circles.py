import geopandas as gpd
from shapely.geometry import Point, Polygon, box, LineString
from shapely.ops import unary_union
import numpy as np
import folium
import math
import os
import time # For adding small delays for UI responsiveness if prints are too fast
from python_tsp.heuristics import solve_tsp_simulated_annealing
from python_tsp.distances import great_circle_distance_matrix
import pulp # For Integer Linear Programming

# --- Configuration & Parameters ---
DRONE_COVERAGE_RADIUS_KM = 10.0 # Drone coverage radius in kilometers

# Helper function to create a circular polygon from a center and radius
def create_circle_geo(center_point: Point, radius_km: float, segments=32) -> Polygon:
    coords = []
    center_lon, center_lat = center_point.x, center_point.y

    # Convert radius from km to degrees for latitude
    radius_deg_lat = radius_km / 111.1 

    # Convert radius from km to degrees for longitude, accounting for latitude
    # This helps make the shape appear more circular on a Web Mercator map.
    # cos(radians(90)) is 0, cos(radians(-90)) is 0. Handle poles.
    if center_lat >= 89.99: # Near North Pole
        # At poles, longitude lines converge; a circle becomes difficult to define this way.
        # Fallback to using latitude degrees for longitude, or handle as a special case (e.g. cap).
        # For simplicity here, we make it behave like at a very high latitude but not exactly the pole.
        effective_lat_for_cos = 89.0 
    elif center_lat <= -89.99: # Near South Pole
        effective_lat_for_cos = -89.0
    else:
        effective_lat_for_cos = center_lat
    
    cos_lat = math.cos(math.radians(effective_lat_for_cos))
    if cos_lat < 0.000001: # Avoid division by zero or extremely small numbers if effective_lat is too close to 90/-90
        radius_deg_lon = radius_deg_lat * 10 # Arbitrary large factor to make it wide, or handle better
    else:
        radius_deg_lon = radius_km / (111.1 * cos_lat)

    for i in range(segments + 1):
        angle_rad = (i / segments) * 2 * math.pi
        # Use the adjusted degree offsets
        dx_deg = radius_deg_lon * math.cos(angle_rad) 
        dy_deg = radius_deg_lat * math.sin(angle_rad) 
        coords.append((center_lon + dx_deg, center_lat + dy_deg))
    return Polygon(coords)

def main():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    grid_geojson_path = os.path.join(project_root, "data", "simulation_outputs", "flood_area_grid.geojson")
    output_circles_geojson_path = os.path.join(project_root, "data", "simulation_outputs", "coverage_circles.geojson")
    output_path_geojson_path = os.path.join(project_root, "data", "simulation_outputs", "coverage_path.geojson")
    output_map_html_path = os.path.join(project_root, "data", "simulation_outputs", "coverage_circles_map.html")

    # --- 1. Load Grid Cells ---
    try:
        grid_gdf = gpd.read_file(grid_geojson_path)
        if grid_gdf.empty:
            print(f"Grid file {grid_geojson_path} is empty. Exiting.")
            return
        if 'cell_id' not in grid_gdf.columns:
            grid_gdf['cell_id'] = range(len(grid_gdf))
        print(f"Loaded {len(grid_gdf)} grid cells from {grid_geojson_path}")
    except Exception as e:
        print(f"Error loading grid GeoJSON: {e}")
        return

    all_cell_ids = grid_gdf['cell_id'].tolist()
    
    # --- 2. Generate All Potential Candidate Circles and Determine Coverage ---
    print("Generating all potential candidate circles and their coverage...")
    candidate_circles_info = [] # List of tuples: (center_point, polygon_shape, set_of_covered_cell_ids)
    
    # For each grid cell, consider its centroid as a potential center for a circle
    for idx, cell_row in grid_gdf.iterrows():
        potential_center = cell_row.geometry.centroid
        potential_circle_poly = create_circle_geo(potential_center, DRONE_COVERAGE_RADIUS_KM)
        
        covered_cell_ids_by_this_candidate = set()
        for cell_id_to_check, check_cell_row in grid_gdf.iterrows():
            if potential_circle_poly.contains(check_cell_row.geometry):
                covered_cell_ids_by_this_candidate.add(check_cell_row['cell_id'])
        
        if covered_cell_ids_by_this_candidate: # Only consider circles that cover at least one cell
            candidate_circles_info.append({
                "center": potential_center,
                "polygon": potential_circle_poly,
                "covers": covered_cell_ids_by_this_candidate,
                "id": f"candidate_circle_{idx}" # Unique ID for PuLP variable
            })
    
    if not candidate_circles_info:
        print("No potential candidate circles could be generated or none cover any cells. Exiting.")
        return
    print(f"Generated {len(candidate_circles_info)} potential candidate circles.")

    # --- 3. Solve Set Cover Problem using Integer Linear Programming (PuLP) ---
    print("Setting up and solving the Set Cover problem using PuLP...")
    start_time_pulp = time.time()

    # Create the LP problem
    prob = pulp.LpProblem("DroneSetCover", pulp.LpMinimize)

    # Decision Variables: One for each candidate circle (binary: 0 if not chosen, 1 if chosen)
    circle_vars = pulp.LpVariable.dicts(
        "Circle", 
        [info["id"] for info in candidate_circles_info], 
        cat=pulp.LpBinary
    )

    # Objective Function: Minimize the sum of chosen circles
    prob += pulp.lpSum([circle_vars[info["id"]] for info in candidate_circles_info]), "TotalCircles"

    # Constraints: Each grid cell must be covered by at least one chosen circle
    for cell_id in all_cell_ids:
        covering_circles_for_this_cell = []
        for info in candidate_circles_info:
            if cell_id in info["covers"]:
                covering_circles_for_this_cell.append(circle_vars[info["id"]])
        
        if covering_circles_for_this_cell: # Only add constraint if the cell can be covered
             prob += pulp.lpSum(covering_circles_for_this_cell) >= 1, f"Cell_{cell_id}_Covered"
        else:
            print(f"Warning: Cell ID {cell_id} cannot be covered by any potential circle. Check grid or radius.")


    # Solve the problem (PuLP will call a solver, CBC is default and usually included)
    # You can specify a solver if needed, e.g., prob.solve(pulp.PULP_CBC_CMD(msg=0)) for no solver messages
    prob.solve() 
    # prob.solve(pulp.PULP_CBC_CMD(msg=True, threads=4)) # Example with CBC options

    pulp_duration = time.time() - start_time_pulp
    print(f"PuLP solver status: {pulp.LpStatus[prob.status]}. Solved in {pulp_duration:.2f} seconds.")

    if prob.status != pulp.LpStatusOptimal:
        print("Optimal solution not found by PuLP. Exiting.")
        if prob.status == pulp.LpStatusInfeasible:
            print("Problem is infeasible - check if all cells can be covered.")
        elif prob.status == pulp.LpStatusUnbounded:
            print("Problem is unbounded - this shouldn't happen for set cover.")
        return

    chosen_circles_centers = []
    chosen_circles_polygons = []
    
    selected_candidate_ids = []
    for info in candidate_circles_info:
        if pulp.value(circle_vars[info["id"]]) == 1: # If this circle was chosen by the solver
            chosen_circles_centers.append(info["center"])
            chosen_circles_polygons.append(info["polygon"])
            selected_candidate_ids.append(info["id"])

    print(f"\nILP algorithm finished. Optimal number of circles found: {len(chosen_circles_polygons)}.")
    print(f"Selected circle candidate IDs: {selected_candidate_ids}")
    
    # --- 4. Save and Visualize --- 
    if not chosen_circles_polygons:
        print("No coverage circles were selected by the ILP solver.")
        return

    circles_gdf = gpd.GeoDataFrame(geometry=chosen_circles_polygons, crs=grid_gdf.crs)
    circles_gdf['circle_id'] = range(len(circles_gdf))

    output_dir = os.path.dirname(output_circles_geojson_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    try:
        circles_gdf.to_file(output_circles_geojson_path, driver="GeoJSON")
        print(f"Successfully saved coverage circles to {output_circles_geojson_path}")
    except Exception as e:
        print(f"Error saving coverage circles GeoJSON: {e}")

    # --- Create and Save Closed Path from Circle Centers (TSP optimization) ---
    path_gdf = None
    if len(chosen_circles_centers) > 1:
        coords_array = np.array([[pt.y, pt.x] for pt in chosen_circles_centers])
        distance_matrix = great_circle_distance_matrix(coords_array)
        permutation, _ = solve_tsp_simulated_annealing(distance_matrix)
        optimized_centers = [chosen_circles_centers[i] for i in permutation]
        
        line_coords = [(pt.x, pt.y) for pt in optimized_centers]
        if len(line_coords) > 1: 
            line_coords.append(line_coords[0]) # Close the loop
            path_line = LineString(line_coords)
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

    # Create Folium map
    if grid_gdf.empty:
        map_center = [0,0] # Fallback
    else:
        map_center_geom = grid_gdf.unary_union.centroid
        map_center = [map_center_geom.y, map_center_geom.x]
    
    m = folium.Map(location=map_center, zoom_start=10, tiles="OpenStreetMap")

    folium.GeoJson(grid_gdf, name='AOI Grid Cells', 
                   style_function=lambda x: {'fillColor': 'gray', 'color': 'black', 'weight':0.5, 'fillOpacity':0.1}
                  ).add_to(m)
    
    folium.GeoJson(circles_gdf, name='Coverage Circles',
                   style_function=lambda x: {'fillColor': 'red', 'color': 'red', 'weight':1, 'fillOpacity':0.3},
                   tooltip=folium.features.GeoJsonTooltip(fields=['circle_id'], aliases=['Circle ID:'])
                  ).add_to(m)

    if path_gdf is not None and not path_gdf.empty:
        folium.GeoJson(path_gdf, name='Coverage Path',
                       style_function=lambda x: {'color': 'blue', 'weight': 3, 'opacity': 0.7}
                      ).add_to(m)

    folium.LayerControl().add_to(m)

    try:
        m.save(output_map_html_path)
        print(f"Successfully saved coverage map to {output_map_html_path}")
    except Exception as e:
        print(f"Error saving map HTML: {e}")

if __name__ == "__main__":
    main() 