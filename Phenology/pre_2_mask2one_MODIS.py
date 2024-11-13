import os
import rasterio
import numpy as np

# List of mask file paths
maskfile_path = r'D:\SetoLab\Phenology\mask'
mask_files = [file for file in os.listdir(maskfile_path) if file.endswith('MODIS.tif')]

# Initialize a variable to store the combined mask
combined_mask = None

# Loop through each file and update the combined mask
# Loop through each file and update the combined mask
for mask_file in mask_files:
    with rasterio.open(os.path.join(maskfile_path, mask_file)) as src:
        # Read the mask data from the file
        mask_data = src.read(1)

        # Initialize the combined mask with the first file's data
        if combined_mask is None:
            combined_mask = mask_data
            combined_meta = src.meta
        else:
            # Combine the masks: pixel is 0 if it's 0 in any of the masks
            combined_mask = np.maximum(combined_mask, mask_data)


# Save the combined mask as a new GeoTIFF file
combined_mask_file = os.path.join(maskfile_path, 'combined_mask_MODIS.tif')
with rasterio.open(combined_mask_file, 'w', **combined_meta) as dest:
    dest.write(combined_mask, 1)
