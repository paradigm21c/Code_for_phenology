# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 23:10:23 2024

@author: parad
"""

import os
from datetime import datetime, timedelta
import rasterio
import numpy as np
import pickle
from scipy.signal import savgol_filter
# from scipy.interpolate import UnivariateSpline
# from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import PchipInterpolator
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor, as_completed
import concurrent.futures
from tqdm import tqdm

def create_complete_date_series(year):
    """
    Create a complete series of daily dates for a given year.
    """
    start_date = datetime(year, 1, 1)
    end_date = datetime(year, 12, 31)
    total_days = (end_date - start_date).days + 1
    return [start_date + timedelta(days=day) for day in range(total_days)]

def insert_nans_for_gaps(pixel_series, actual_dates, complete_dates):
    """
    Insert NaNs in the pixel_series where there are gaps in the actual_dates.
    """
    pixel_series_with_gaps = []
    actual_date_index = 0

    for date in complete_dates:
        if actual_date_index < len(actual_dates) and date == actual_dates[actual_date_index]:
            pixel_series_with_gaps.append(pixel_series[actual_date_index])
            actual_date_index += 1
        else:
            pixel_series_with_gaps.append(np.nan)

    return np.array(pixel_series_with_gaps)

def fill_gaps(pixel_series_with_gaps):
    """
    Fill gaps in the pixel_series using 'Piecewise Cubic Hermite Interpolating Polynomial (PCHIP)'.
        Method: Ensures monotonicity between data points, preventing overshoots in the interpolated curve.
        Use Case: Best for data where maintaining the shape of the curve (like increasing or decreasing trends) is important.
    Assumes NaNs represent gaps.
    """
    valid_mask = ~np.isnan(pixel_series_with_gaps)
    x_valid = np.where(valid_mask)[0]
    y_valid = pixel_series_with_gaps[valid_mask]

    if len(x_valid) < 2:
        return pixel_series_with_gaps

    try:
        spline = PchipInterpolator(x_valid, y_valid)
        return spline(np.arange(len(pixel_series_with_gaps)))
    except Exception as e:
        print(f"Error during interpolation: {e}")
        # Handle the error or return the original series with gaps
        return pixel_series_with_gaps

def process_pixel(pixel_series, temp_date):
    year_list = np.unique([date.year for date in temp_date])
    complete_dates = create_complete_date_series(year_list[0])
    
    pixel_series_with_gaps = insert_nans_for_gaps(pixel_series, temp_date, complete_dates)
    pixel_series_with_nogaps=fill_gaps(pixel_series_with_gaps)             
                        
    # Updated function to process a single pixel's time series
    try:
        # Identify the 'background' EVI2 value (10th percentile of snow-free EVI2 values)
        background_evi2 = np.percentile(pixel_series_with_nogaps, 10)

        # Replace EVI2 values below the background value with the background EVI2
        pixel_series = np.where(pixel_series_with_nogaps < background_evi2, background_evi2, pixel_series_with_nogaps)

        # Smooth time series
        smoothed_evi2 = savgol_filter(pixel_series, window_length=8, polyorder=2)
        smoothed_evi2 = savgol_filter(smoothed_evi2, window_length=8, polyorder=6)
        
        # Calculate min, max, dates for each year in the pixel series
        return calculate_min_max_dates(smoothed_evi2, complete_dates)
    except Exception as e:
        print(f"Error processing pixel: {e}")
        return None 

def calculate_evi2(nir, red):
    evi2 = 2.5 * (nir - red) / (nir + 2.4 * red + 1)
    # Set evi2 to nan where it is out of the [0, 1] range
    evi2[(evi2 < 0) | (evi2 > 1)] = np.nan
    return evi2

def extract_date_from_filename(filename):
    # Extracting date assuming format 'merge_YYYYMMDD_PF.tif'
    date_str = os.path.basename(filename).split('_')[1]
    return datetime.strptime(date_str, '%Y%m%d')

def calculate_min_max_dates(smoothed_evi2, complete_dates):
    dates = complete_dates
    years = np.array([date.year for date in dates])
    year=years[0]
    results = []

    year_series = smoothed_evi2
    year_dates = np.array(dates)

    # Find indices of min and max values
    min_index = np.argmin(year_series)
    max_index = np.argmax(year_series)

    # Corresponding min and max values
    min_val = year_series[min_index]
    max_val = year_series[max_index]
    
    # Ensure criteria for valid growth cycle are met
    evi2_diff = max_val - min_val
    if evi2_diff < 0.1:
        year_result = {
            'year': year,
            'min_value': np.nan,
            'min_date': np.nan,
            'max_value': np.nan,
            'max_date': np.nan,
            'val_15_value': np.nan,
            'greenup_date': np.nan,
            'dormancy_date': np.nan,
            'val_50_value': np.nan,
            'midgreenup_date': np.nan,
            'midgreendown_date': np.nan,
            'val_90_value': np.nan,
            'maturity_date': np.nan,
            'senescence_date': np.nan,
            'val_sos': np.nan,
            'sos_date' : np.nan,
            'val_eos' : np.nan,
            'eos_date' : np.nan,
        }
        results.append(year_result)
    else:
        # Corresponding dates for min and max values
        min_date = year_dates[min_index]
        max_date = year_dates[max_index]
        
        # Process left and right parts of the series
        year_dates_left, year_series_left = year_dates[:max_index+1], year_series[:max_index+1]
        year_dates_right, year_series_right = year_dates[max_index:], year_series[max_index:]
        
        # Find indices of min left and right values
        min_index_left = np.argmin(year_series_left)
        min_index_right = np.argmax(year_series_right)
    
        # Corresponding min left and right values
        min_val_left = year_series_left[min_index_left]
        min_val_right = year_series_right[min_index_right]
    
        # Calculate threshold and find the first crossing date
        amp = max_val - min_val
        amp_left = max_val - min_val_left
        amp_right = max_val - min_val_right
        
        val_15 = min_val + amp * 0.15
        val_50 = min_val + amp * 0.50
        val_90 = min_val + amp * 0.90
        
        val_sos = min_val_left + amp_left * 0.20
        val_eos = min_val_right + amp_right * 0.20 
        
        val_15_dates = year_dates[year_series >= val_15]
        val_50_dates = year_dates[year_series >= val_50]
        val_90_dates = year_dates[year_series >= val_90]
        
        Sos_dates = year_dates_left[year_series_left >= val_sos]
        Eos_dates = year_dates_right[year_series_right >= val_eos]
        
        Greenup_date = val_15_dates[0] if len(val_15_dates) > 0 else None
        MidGreenup_date = val_50_dates[0] if len(val_50_dates) > 0 else None
        Maturity_date = val_90_dates[0] if len(val_90_dates) > 0 else None
        
        Dormancy_date = val_15_dates[-1] if len(val_15_dates) > 0 else None
        MidGreendown_date = val_50_dates[-1] if len(val_50_dates) > 0 else None
        Senescence_date = val_90_dates[-1] if len(val_90_dates) > 0 else None
        
        Sos_date = Sos_dates[0] if len(Sos_dates) > 0 else None
        Eos_date = Eos_dates[0] if len(Eos_dates) > 0 else None
    
    
        # Store the results for this year
        year_result = {
            'year': year,
            'min_value': min_val,
            'min_date': min_date,
            'max_value': max_val,
            'max_date': max_date,
            
            'val_15_value': val_15,
            'greenup_date': Greenup_date,
            'dormancy_date': Dormancy_date,
            
            'val_50_value': val_50,
            'midgreenup_date': MidGreenup_date,
            'midgreendown_date': MidGreendown_date,
            
            'val_90_value': val_90,
            'maturity_date': Maturity_date,
            'senescence_date': Senescence_date,
            
            'val_sos': val_sos,
            'sos_date' : Sos_date,
            'val_eos' : val_eos,
            'eos_date' : Eos_date,
            
        }
        results.append(year_result)

    return results

def process_evi2(file_path, mask):
    with rasterio.open(file_path) as src:
        nir = src.read(4)
        red = src.read(3)
        nir[~mask] = 0
        red[~mask] = 0
        
        # Band 4 is NIR and band2 3 is Red
        nir = nir.astype(float) / 10000
        red = red.astype(float) / 10000
        
        # Convert values outside the range [0, 1] to NaN
        nir = np.where((nir >= 0) & (nir <= 1), nir, np.nan)
        red = np.where((red >= 0) & (red <= 1), red, np.nan)

        return calculate_evi2(nir, red)

def process_files(target_directory, temp_file_list, mask, temp_date, product, year, data_save_path):
    with rasterio.open(os.path.join(target_directory, temp_file_list[0])) as src:
        height, width = src.shape
    # Initialize an empty list to hold the EVI2 maps
    evi2_maps = [None] * len(temp_file_list)

    # Use ThreadPoolExecutor to parallelize file processing
    with ThreadPoolExecutor() as executor:
        # Create a future for each file
        future_to_index = {executor.submit(process_evi2, os.path.join(target_directory, file), mask): i for i, file in enumerate(temp_file_list)}

        for future in as_completed(future_to_index):
            index = future_to_index[future]
            try:
                evi2_maps[index] = future.result()
            except Exception as exc:
                print(f'File {temp_file_list[index]} generated an exception: {exc}')

    # Stack the processed EVI2 maps
    stacked_evi2_maps = np.stack(evi2_maps)
   
    processed_data = {}
    futures = {}

    # Calculate total number of pixels to process
    total_pixels = np.sum(mask)  

    # Initialize ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=10) as executor:
        for x in range(height):
            for y in range(width):
                if mask[x, y]:
                    pixel_series = stacked_evi2_maps[:, x, y]                  
                    future = executor.submit(process_pixel, pixel_series, temp_date)
                    futures[future] = (x, y)

        # Using tqdm to show progress bar
        for future in tqdm(concurrent.futures.as_completed(futures), total=total_pixels, desc="Processing Pixels"):
            x, y = futures[future]
            result = future.result()
            if result is not None:
                processed_data[(x, y)] = result

    # Save processed_data to a file
    output_file_path = os.path.join(data_save_path, f'processed_data_{product}_{year}.pkl')
    with open(output_file_path, 'wb') as f:
        pickle.dump(processed_data, f)

    print(f"Data saved to {output_file_path}")
    
def main():
    products =['SDC'] #'MODIS']#,'cuFSDAF'] #'Landsat']
    
    for product in products:
        if product == 'MODIS':
            target_directory = r'D:\SetoLab\MODIS\MCD43A4\Pre_Fusion_Smoothed'
        elif product == 'Landsat':
            target_directory = r'D:\SetoLab\One30m\Pre_Fusion_500m'
        elif product == 'cuFSDAF':    
            target_directory = r'D:\SetoLab\MODIS\MCD43A4\Fusion_500m'
        elif product == 'SDC':    
            target_directory = r'C:\Temp\MODIS\SDC500\Fusion_500m'             
    
        data_save_path = r'D:\SetoLab\Phenology\data_cal'
        masks_path = r'D:\SetoLab\Phenology\mask'
        
        file_list = [file for file in os.listdir(target_directory) if file.endswith('.tif')]
        mask_path_file = os.path.join(masks_path, f'combined_mask_{product}.tif')  # Adjusted to use temp_name

        dates = [extract_date_from_filename(f) for f in file_list]
        years = np.array([date.year for date in dates])
        
        for year in [2016,2017,2018,2019,2020,2021,2022]: 
            # Inside your loop
            temp_file_list = [file for file, y in zip(file_list, years) if y == year]
            temp_date = [extract_date_from_filename(f) for f in temp_file_list]
            
            with rasterio.open(mask_path_file) as src:
                mask = src.read(1).astype(bool)
            
            process_files(target_directory, temp_file_list, mask, temp_date, product, year, data_save_path)


if __name__ == '__main__':
    main()