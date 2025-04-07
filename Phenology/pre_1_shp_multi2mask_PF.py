import os
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from concurrent.futures import ThreadPoolExecutor, as_completed

# Paths
maskfile_path = 'D:/SetoLab/Phenology/mask/parks_arc_PF'
shapefile_path = 'D:/SetoLab/Phenology/QGIS/ArcGIS_map'
geotiff_path = 'C:/Temp/PF/merged_20230801_PF.tif'  # 6-band Landsat TIFF file

# Load the shapefile
shapefile_full_path = os.path.join(shapefile_path, "US_park_UTM_18N_JW_NYC.shp")
shapes = gpd.read_file(shapefile_full_path)

# Initialize a counter for polygons without names
nan_counter = 1

# Process each polygon in the shapefile
for index, polygon in shapes.iterrows():
    # Skip if the polygon's geometry is None
    if polygon.geometry is None:
        print(f"Skipping polygon at index {index} due to missing geometry.")
        continue
    
    # Ensure there's a 'name' column in your shapefile. Replace 'name' below with the correct column name if different.
    if 'NAME' in polygon and isinstance(polygon['NAME'], str):
        # Replace spaces with underscores in the polygon name
        polygon_name = polygon['NAME'].replace(' ', '_')
        polygon_name = polygon_name.replace("'", '')
        polygon_name = polygon_name.replace(".", '')
        polygon_name = polygon_name.replace("-", '')
    # else:
        # print(f"No name found for polygon at index {index}, skipping file save.")
        # continue  # Skip this polygon if no name is provided
        
    else:
        # Assign a default name like nan1, nan2, etc.
        polygon_name = f"nan{nan_counter}"
        nan_counter += 1
    
    with rasterio.open(geotiff_path) as src:
        # Mask operation for the first band using the current polygon
        out_image, out_transform = mask(src, [polygon.geometry], all_touched=False, crop=False, indexes=1)
        out_meta = src.meta

        # Update metadata for masked output
        out_meta.update({
            "driver": "GTiff",
            "height": src.height,
            "width": src.width,
            "transform": out_transform,
            "count": 1,
            "nodata": 9999,
        })

        # Check if there are any masked pixels
        if out_image.sum() > 0:
            # Define output filename using the sanitized polygon name
            output_filename = f"{polygon_name}_mask.tif"
            output_path = os.path.join(maskfile_path, output_filename)

            # Write the mask file
            with rasterio.open(output_path, 'w', **out_meta) as dest:
                dest.write(out_image, 1)
        else:
            print(f"No masked pixels for polygon named {polygon_name}, skipping file save.")

print("Finished creating masks for all polygons.")

