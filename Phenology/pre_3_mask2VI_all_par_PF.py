import os
from datetime import datetime
import rasterio
import numpy as np
import pickle
from scipy.signal import savgol_filter
# from scipy.interpolate import UnivariateSpline
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor, as_completed
import concurrent.futures
from tqdm import tqdm

'''       
def process_pixel(pixel_series, file_list):
    try:
        # Identify the 'background' EVI2 value (10th percentile of snow-free EVI2 values)
        background_evi2 = np.percentile(pixel_series, 10)

        # Replace EVI2 values below the background value with the background EVI2
        pixel_series = np.where(pixel_series < background_evi2, background_evi2, pixel_series)

        # Smooth time series using penalized cubic smoothing spline
        x = np.arange(len(pixel_series))  # Assuming equal spacing
        spline = UnivariateSpline(x, pixel_series, k=3, s=some_smoothing_factor)  # Adjust 's' as needed
        smoothed_evi2 = spline(x)

        # Calculate min, max, dates for each year in the pixel series
        return calculate_min_max_dates(smoothed_evi2, file_list)
    except Exception as e:
        print(f"Error processing pixel: {e}")
        return None
'''

def process_pixel(pixel_series, file_list):
    # Updated function to process a single pixel's time series
    try:
        # Identify the 'background' EVI2 value (10th percentile of snow-free EVI2 values)
        background_evi2 = np.percentile(pixel_series, 10)

        # Replace EVI2 values below the background value with the background EVI2
        pixel_series = np.where(pixel_series < background_evi2, background_evi2, pixel_series)

        # Smooth time series
        smoothed_evi2 = savgol_filter(pixel_series, window_length=8, polyorder=2)
        smoothed_evi2 = savgol_filter(smoothed_evi2, window_length=8, polyorder=6)
        
        # Calculate min, max, dates for each year in the pixel series
        return calculate_min_max_dates(smoothed_evi2, file_list)
    except Exception as e:
        print(f"Error processing pixel: {e}")
        return None 

def calculate_evi2(nir, red):
    return 2.5 * (nir - red) / (nir + 2.4 * red + 1)

def extract_date_from_filename(filename):
    # Extracting date assuming format 'merge_YYYYMMDD_PF.tif'
    date_str = os.path.basename(filename).split('_')[1]
    return datetime.strptime(date_str, '%Y%m%d')

def calculate_min_max_dates(pixel_series, file_list):
    dates = [extract_date_from_filename(f) for f in file_list]
    years = np.array([date.year for date in dates])

    results = []

    for year in np.unique(years):
        year_mask = years == year
        year_series = pixel_series[year_mask]
        year_dates = np.array(dates)[year_mask]
        
        if len(year_series) == 0:
            continue

        # Find indices of min and max values
        min_index = np.argmin(year_series)
        max_index = np.argmax(year_series)

        # Corresponding min and max values
        min_val = year_series[min_index]
        max_val = year_series[max_index]
        
        # Ensure criteria for valid growth cycle are met
        evi2_diff = max_val - min_val
        if evi2_diff < 0.1 or evi2_diff < 0.35 * np.ptp(pixel_series):
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
            }
            results.append(year_result)
            continue

        # Corresponding dates for min and max values
        min_date = year_dates[min_index]
        max_date = year_dates[max_index]

        # Calculate threshold and find the first crossing date
        amp = max_val - min_val
        val_15 = min_val + amp * 0.15
        val_50 = min_val + amp * 0.50
        val_90 = min_val + amp * 0.90
        val_15_dates = year_dates[year_series >= val_15]
        val_50_dates = year_dates[year_series >= val_50]
        val_90_dates = year_dates[year_series >= val_90]
        
        Greenup_date = val_15_dates[0] if len(val_15_dates) > 0 else None
        MidGreenup_date = val_50_dates[0] if len(val_50_dates) > 0 else None
        Maturity_date = val_90_dates[0] if len(val_90_dates) > 0 else None
        
        Dormancy_date = val_15_dates[-1] if len(val_15_dates) > 0 else None
        MidGreendown_date = val_50_dates[-1] if len(val_50_dates) > 0 else None
        Senescence_date = val_90_dates[-1] if len(val_90_dates) > 0 else None


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
        }
        results.append(year_result)

    return results

def process_evi2(file_path, mask):
    with rasterio.open(file_path) as src:
        nir = src.read(4)
        red = src.read(3)
        nir[~mask] = 0
        red[~mask] = 0
        return calculate_evi2(nir*0.0001, red*0.0001)

def process_files(files_path, file_list, mask):
    with rasterio.open(os.path.join(files_path, file_list[0])) as src:
        height, width = src.shape
    # Initialize an empty list to hold the EVI2 maps
    evi2_maps = [None] * len(file_list)

    # Use ThreadPoolExecutor to parallelize file processing
    with ThreadPoolExecutor() as executor:
        # Create a future for each file
        future_to_index = {executor.submit(process_evi2, os.path.join(files_path, file), mask): i for i, file in enumerate(file_list)}

        for future in as_completed(future_to_index):
            index = future_to_index[future]
            try:
                evi2_maps[index] = future.result()
            except Exception as exc:
                print(f'File {file_list[index]} generated an exception: {exc}')

    # Stack the processed EVI2 maps
    stacked_evi2_maps = np.stack(evi2_maps)

    processed_data = {}
    futures = {}

    # Calculate total number of pixels to process
    total_pixels = np.sum(mask)

    # Initialize ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=30) as executor:
        for x in range(height):
            for y in range(width):
                if mask[x, y]:
                    pixel_series = stacked_evi2_maps[:, x, y]
                    future = executor.submit(process_pixel, pixel_series, file_list)
                    futures[future] = (x, y)

        # Using tqdm to show progress bar
        for future in tqdm(concurrent.futures.as_completed(futures), total=total_pixels, desc="Processing Pixels"):
            x, y = futures[future]
            result = future.result()
            if result is not None:
                processed_data[(x, y)] = result

    return processed_data
    
def main():
    products =['cuFSDAF'] #'Landsat']
    
    for product in products:
        if product == 'MODIS':
            target_directory = r'D:\SetoLab\MODIS\MCD43A4\Pre_Fusion_Smoothed'
        elif product == 'Landsat':
            target_directory = r'D:\SetoLab\One30m\Pre_Fusion_500m'
        elif product == 'cuFSDAF':    
            target_directory = r'D:\SetoLab\MODIS\MCD43A4\Fusion_500m'    
    
    data_save_path = r'D:\SetoLab\Phenology\data_cal'
    masks_path = r'D:\SetoLab\Phenology\mask'
    
    file_list = [file for file in os.listdir(target_directory) if file.endswith('.tif')]
    mask_path_file = os.path.join(masks_path, f'combined_mask_{product}.tif')  # Adjusted to use temp_name

    dates = [extract_date_from_filename(f) for f in file_list]
    years = np.array([date.year for date in dates])
    
    for year in np.unique(years):
        # Inside your loop
        temp_file_list = [file for file, y in zip(file_list, years) if y == year]

        with rasterio.open(mask_path_file) as src:
            mask = src.read(1).astype(bool)
        
        processed_data = process_files(target_directory, temp_file_list, mask)

        # Save processed_data to a file
        output_file_path = os.path.join(data_save_path, f'processed_data_{product}_{year}.pkl')
        with open(output_file_path, 'wb') as f:
            pickle.dump(processed_data, f)
        print(f"Data saved to {output_file_path}")

if __name__ == '__main__':
    main()