# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 22:14:36 2023

@author: parad
"""

import os
from collections import defaultdict
from datetime import datetime, timedelta
from multiprocessing import Pool, cpu_count
import rasterio
from rasterio.warp import reproject, calculate_default_transform
from rasterio.enums import Resampling
from rasterio.crs import CRS
from osgeo import gdal
from concurrent.futures import ProcessPoolExecutor
import shutil

gdal.UseExceptions()

def delete_folder(folder_path):
    try:
        # Delete the folder and its contents
        shutil.rmtree(folder_path)
        print(f"Deleted folder: {folder_path}")
    except Exception as e:
        print(f"Error deleting folder: {str(e)}")

def stack_bands(date_files_pair, output_dir):
    date, file_list = date_files_pair
    
    # Define the desired order of bands
    desired_order = [3, 4, 1, 2, 5, 6] 

    # Reorder the file_list based on the desired order of bands
    # Extract band numbers from the filenames and match them to the desired order for sorting
    file_list_sorted = sorted(file_list, key=lambda x: desired_order.index(int(x.split('_')[-1].split('.')[0][1:])))

    # Adjust sorting logic for the specific filename format
    file_list_sorted = sorted(file_list, key=lambda x: desired_order.index(int(x.split('_b')[-1].split('.')[0])))

    # Read the metadata from the first file
    with rasterio.open(file_list_sorted[0]) as src0:
        meta = src0.meta

    # Update meta to reflect the number of layers and datatype
    meta.update(count=len(file_list_sorted), dtype='int16')

    # Write the stacked bands to a new file in the desired order
    output_path = os.path.join(output_dir, f'stacked_{date}_MODIS.tif')
    with rasterio.open(output_path, 'w', **meta) as dst:
        for id, layer in enumerate(file_list_sorted, start=1):
            with rasterio.open(layer) as src1:
                dst.write_band(id, src1.read(1).astype('int16'))

    return f"Stacking for {date} completed."

def stack_helper(args):
    return stack_bands(*args)


def doy_to_date(year: int, doy: int) -> str:
    """Convert a given year and day-of-year to YYYYMMDD format."""
    date = datetime(year, 1, 1) + timedelta(doy - 1)
    return date.strftime('%Y%m%d')
    

def reproject_MODIS(args):
    """
    Reprojects and resamples a raster dataset.
    
    Parameters:
    - input_path: path to the input raster file
    - output_path: path where the reprojected and resampled raster will be saved
    - new_crs: the new coordinate reference system (can be an EPSG code or a proj4 string)
    - new_resolution: the new resolution for the raster
    """
    
    input_path, output_path, new_crs = args
       
    # Define the CRS for MODIS Sinusoidal (https://spatialreference.org/ref/sr-org/modis-sinusoidal-3/)
    src_crs = 'PROJCS["MODIS Sinusoidal",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Sinusoidal"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",0.0],PARAMETER["semi_major",6371007.181],PARAMETER["semi_minor",6371007.181],UNIT["m",1.0],AUTHORITY["SR-ORG","6974"]]'

    with rasterio.open(input_path) as src:
        # Use the CRS from the source GeoTIFF file
        #src_crs = src_crs
        
        transform, width, height = calculate_default_transform(
            src_crs, new_crs, src.width, src.height, *src.bounds)
    
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': new_crs,
            'transform': transform,
            'width': width,
            'height': height
        })

        with rasterio.open(output_path, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src_crs,
                    src_nodata=None,
                    dst_transform=transform,                   
                    dst_crs=new_crs,
                    dst_nodata = None,
                    resampling=Resampling.cubic_spline
                )

def resample_MODIS(args):
    input_path, output_path, new_resolution = args
    x_res = new_resolution[0]
    y_res = new_resolution[1]

    gdal.Warp(output_path, input_path, xRes=x_res, yRes=y_res, resampleAlg=gdal.GRA_CubicSpline)
                   

def extract_wrapper(hdf_files_by_date, band_names_to_extract, extracted_paths):
    extract_bands_from_hdf(hdf_files_by_date, band_names_to_extract, extracted_paths)
    
def reproject_data(input_directory, output_directory, new_crs):
    files_to_process = [(os.path.join(input_directory, file),
                         os.path.join(output_directory, f"{file.replace('extracted', 'reprojected')}"),
                         new_crs)
                        for file in os.listdir(input_directory) if file.endswith('.tif')]

    num_cpus = 20 #cpu_count()
    
    with Pool(processes=num_cpus) as pool:
        pool.map(reproject_MODIS, files_to_process)
        
def resample_data(input_directory, output_directory, new_resolution):
    files_to_process = [(os.path.join(input_directory, file),
                         os.path.join(output_directory, f"{file.replace('stacked', 'resampled')}"),
                         new_resolution)
                        for file in os.listdir(input_directory) if file.endswith('.tif')]

    num_cpus = 24 #cpu_count()
    
    with Pool(processes=num_cpus) as pool:
        pool.map(resample_MODIS, files_to_process)
    
def main():
     # Step 1: Stacking TIFF
    base_path_directory = r'F:/SDC/dd'
    stacked_directory = r'C:/Temp/MODIS/Stacked'
    # stacked_directory = r'F:/SDC_Stacked/temp'

    if not os.path.exists(stacked_directory):
        os.makedirs(stacked_directory)

    files_by_date = defaultdict(list)
    
    for year_dir in os.listdir(base_path_directory):
        year_path = os.path.join(base_path_directory, year_dir)
        if os.path.isdir(year_path):
            for filename in os.listdir(year_path):
                if filename.endswith('.tif'):
                    parts = filename.split('_')
                    year_doy = parts[3]
                    year, doy = int(year_doy[:4]), int(year_doy[4:])
                    date = doy_to_date(year, doy)
                    files_by_date[date].append(os.path.join(year_path, filename))
    
    with ProcessPoolExecutor(max_workers=10) as executor:
        tasks = [(date_files, stacked_directory) for date_files in files_by_date.items()]
        results = executor.map(stack_helper, tasks)

    for r in results:
        print(r)

    print("All band stacking completed.")
    
    # Step 2: Reproject
    reprojected_directory = r"C:/Temp/MODIS/Reprojected"
    new_crs = CRS.from_epsg(32617)  # Example: UTM Zone 18N
    
    if not os.path.exists(reprojected_directory):
        os.makedirs(reprojected_directory)
    
    reproject_data(stacked_directory, reprojected_directory, new_crs)

    print("Reprojection completed.")
    delete_folder(stacked_directory)
    
     # Step 3: Resample
    resampled_directory = r"C:/Temp/MODIS/Resampled"
    new_resolution = (480, 480)  # Example: 30m x 30m
    
    if not os.path.exists(resampled_directory):
        os.makedirs(resampled_directory)
    
    resample_data(reprojected_directory, resampled_directory, new_resolution)

    print("Resampling completed.")
    # delete_folder(reprojected_directory)

if __name__ == '__main__':
    main()