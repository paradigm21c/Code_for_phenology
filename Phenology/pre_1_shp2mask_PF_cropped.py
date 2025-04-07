import os
import geopandas as gpd
import rasterio
from rasterio.mask import mask

# Save mask path
maskfile_path = r'D:\SetoLab\Phenology\mask'

# Load the shapefile
shapefile_path = r'D:\SetoLab\NAC_ECM\SHP'
shp_files = [file for file in os.listdir(shapefile_path) if file.endswith('.shp')]

reference_directory = r"D:\SetoLab\code\PF\crop"
reference_files = [os.path.join(reference_directory, f) for f in os.listdir(reference_directory) if f.endswith('.tif')]

for shp_file in shp_files:
    shapes = gpd.read_file(os.path.join(shapefile_path,shp_file))
    
    # Load the GeoTIFF file
    for reference_file in reference_files:
        save_name = os.path.basename(reference_file)
        
        with rasterio.open(reference_file) as src:
            # Select the band you want to mask (e.g., band 1)
            selected_band = src.read(1)
        
            # Create a mask for each geometry in the shapefile
            out_image, out_transform = mask(src, shapes.geometry, crop=False, indexes=1)
            out_meta = src.meta
            out_image[out_image !=0]=1
        
        # Update metadata for the mask file
        out_meta.update({"driver": "GTiff",
                         "count": 1,  # Only one band
                         "height": src.shape[0],
                         "width": src.shape[1],
                         "transform": out_transform})
        
        output_name = shp_file.replace('.shp',f'_{save_name}')
        # Write the mask file
        with rasterio.open(os.path.join(maskfile_path,f"{output_name }"), 'w', **out_meta) as dest:
            dest.write(out_image, 1)  # Write only the selected band