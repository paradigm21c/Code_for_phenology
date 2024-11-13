import os
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject, calculate_default_transform
from rasterio.crs import CRS
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor

# Bands you want to extract by sensor type
# Landsat 4/5/7 = B_06, QA_pixel
# Landsat 8/9 = B_10, QA_pixel

bands_to_extract_map = {
    'LT04': [5,2],  
    'LT05': [5,2],  
    'LE07': [5,2],  
    'LC08': [5,2],
    'LC09': [5,2],    
    # Add any other sensor type mappings here
    }

def extract_and_reproject(input_path, output_path, bands_map, new_crs, new_resolution):
    with rasterio.open(input_path) as src:
        # Extract specific bands
        sensor_type = os.path.basename(input_path).split('_')[-1].split('.')[0]
        bands = bands_map.get(sensor_type, [])
        band_data = [src.read(b) for b in bands]
        
        # Prepare metadata for reprojection
        transform, width, height = rasterio.warp.calculate_default_transform(
            src.crs, new_crs, src.width, src.height, *src.bounds, resolution=new_resolution)
        meta = src.meta.copy()
        meta.update({
            'crs': new_crs,
            'transform': transform,
            'width': width,
            'height': height,
            'count': len(bands),
            'nodata': 0
        })

        # Save the reprojected & extracted bands
        with rasterio.open(output_path, 'w', **meta) as dst:
            for idx, data in enumerate(band_data, 1):
                if idx == len(band_data):  # if it's the last band
                    resampling_method = Resampling.nearest
                else:
                    resampling_method = Resampling.bilinear
                    
                reproject(
                    source=data, #*0.00341802 + 149.0)*100,
                    destination=rasterio.band(dst, idx),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    src_nodata=None,
                    dst_transform=transform,
                    dst_crs=new_crs,
                    dst_nodata = 0,
                    resampling=resampling_method
                )

    print(f"Processed {input_path} to {output_path}")

def main():
    input_directory = r"E:\Jay\Satellite\Landsat_ARD\Stacked_ST" # "C:\Temp\ST_NYC"
    output_directory = r"C:\Temp\Reprojected_Extracted_ST"
    # reference_file = "D:\SetoLab\Landsat_ARD\stacked_20130414_L30.tif"
    # src = rasterio.open(reference_file)
    new_crs = CRS.from_epsg(32618)
    new_resolution = (30, 30)

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    args_list = [
        (os.path.join(input_directory, filename), 
         os.path.join(output_directory, filename.replace('stacked_','reprojected_')), 
         bands_to_extract_map, 
         new_crs, 
         new_resolution) 
        for filename in os.listdir(input_directory) if filename.endswith(".tif")
    ]

    num_cpus = cpu_count()

    # Use ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=5) as executor:
        # Instead of using lambda, unpack the tuple arguments directly in the map function
        executor.map(extract_and_reproject, 
                     [args[0] for args in args_list],
                     [args[1] for args in args_list],
                     [args[2] for args in args_list],
                     [args[3] for args in args_list],
                     [args[4] for args in args_list])

    print("Processing completed.")

if __name__ == '__main__':
    main()