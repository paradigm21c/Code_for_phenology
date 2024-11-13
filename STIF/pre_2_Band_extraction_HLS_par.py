import os
import rasterio
from rasterio.crs import CRS
from concurrent.futures import ProcessPoolExecutor
from functools import partial

def extract_bands_from_file(filepath, output_dir, bands_map):
    """Extract specific bands from a tif file and save to a new tif."""
    with rasterio.open(filepath) as src:
        # Extract metadata
        meta = src.meta
        meta['crs']= CRS.from_epsg(32619)

        # Determine the sensor type from the filename
        sensor_type = os.path.basename(filepath).split('_')[-1].split('.')[0]
        
        # Get the bands to extract based on the sensor type
        bands = bands_map.get(sensor_type, [])
        meta['count'] = len(bands)
        
        # Create output file name based on input name and bands extracted
        base_name = os.path.basename(filepath)
        new_filename = f"{base_name.replace('merged','extracted')}"
        output_path = os.path.join(output_dir, new_filename)
        
        with rasterio.open(output_path, 'w', **meta) as dest:
            for idx, band_num in enumerate(bands, 1):
                band_data = src.read(band_num)
                dest.write(band_data, idx)
                
    print(f"Bands {bands} extracted from {filepath} to {output_path}")


def process_file(filepath, output_dir, bands_map):
    extract_bands_from_file(filepath, output_dir, bands_map)

def main():
    # Directory containing the original tif files and directory for output files
    path_directory = r'C:\Temp\HLS\Merged'
    path_output = r'C:\Temp\HLS\ExtractedBands'

    # Bands you want to extract by sensor type
    bands_to_extract_map = {
        'L30': [2, 3, 4, 5, 6, 7, 11],  # Adjust as necessary
        'S30': [2, 3, 4, 13, 11, 12, 14]   # Adjust as necessary
    }    
    
    if not os.path.exists(path_output):
        os.makedirs(path_output)

    # Create a partial function with predefined arguments for process_file
    process_file_partial = partial(process_file, output_dir=path_output, bands_map=bands_to_extract_map)

    # Use ProcessPoolExecutor to parallelize the extraction process
    with ProcessPoolExecutor(max_workers=20) as executor:
        for filename in os.listdir(path_directory):
            if filename.endswith('.tif'):
                filepath = os.path.join(path_directory, filename)
                executor.submit(process_file_partial, filepath)

if __name__ == '__main__':
    main()