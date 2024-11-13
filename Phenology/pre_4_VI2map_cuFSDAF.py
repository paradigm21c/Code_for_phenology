import pickle
import rasterio
import os
import numpy as np
from datetime import datetime


def load_processed_data(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)
        
def date_to_doy(date):
    """Converts a datetime object to Day of Year."""
    # Check if the date is a NaN value (float)
    if isinstance(date, float) and np.isnan(date):
        return np.nan  # Return NaN if the date is NaN
    elif isinstance(date, datetime):  # Check if it's a datetime object
        return date.timetuple().tm_yday
    else:
        raise TypeError("Expected datetime object, got {}".format(type(date).__name__))
  
def create_maps():
    save_name = 'SDC'
    
    path_cal = r'D:\SetoLab\Phenology\data_cal'
    path_map = r'D:\SetoLab\Phenology\map'
    reference_directory = r"D:\SetoLab\Phenology\mask"
    reference_file = f'D:\SetoLab\Phenology\mask\combined_mask_{save_name}.tif'
    
    start_year, end_year = 2000,2022
    
    geotiff_path = os.path.join(reference_directory, reference_file)

    save_path =os.path.join(path_map, save_name)

    if not os.path.exists(save_path):
        os.makedirs(save_path)
        
    with rasterio.open(geotiff_path) as src:
        height = src.height
        width = src.width
        profile = src.profile
    
    # Assuming you know the range of years to process
    for year in range(start_year, end_year + 1):
        processed_data_file = os.path.join(path_cal, f'processed_data_{save_name}_{year}.pkl')
        if not os.path.exists(processed_data_file):
            continue  # Skip if the file doesn't exist for this year

        processed_data = load_processed_data(processed_data_file)

        # Initialize dictionaries for different variables
        maps = {
            'max_value': np.full((height, width), np.nan), 'min_value': np.full((height, width), np.nan),
            'max_date': np.full((height, width), np.nan), 'min_date': np.full((height, width), np.nan), 
            'val_15_value': np.full((height, width), np.nan), 'val_50_value': np.full((height, width), np.nan), 'val_90_value': np.full((height, width), np.nan), 
            'greenup_date': np.full((height, width), np.nan), 'dormancy_date': np.full((height, width), np.nan),
            'midgreenup_date': np.full((height, width), np.nan), 'midgreendown_date': np.full((height, width), np.nan), 
            'maturity_date': np.full((height, width), np.nan), 'senescence_date': np.full((height, width), np.nan),
            'val_sos': np.full((height, width), np.nan),  'sos_date' : np.full((height, width), np.nan),
            'val_eos' : np.full((height, width), np.nan), 'eos_date' : np.full((height, width), np.nan),
        }

        for (x, y), pixel_list in processed_data.items():
            for pixel_data in pixel_list:
                for var in maps:
                    if var in pixel_data:
                        if var.endswith("_date"):
                            try:
                                new_value = date_to_doy(pixel_data[var])
                                maps[var][x, y] = new_value
                            except Exception as e:
                                print(f"Error processing date {pixel_data[var]}: {e}")
                        else:
                            maps[var][x, y] = pixel_data[var]



        for var_name, map_data in maps.items():
            output_filename = os.path.join(save_path, f"EVI2_{save_name}_{var_name}_{year}.tif")
            profile.update(dtype=rasterio.float32)  # Update the profile with the correct data type
            with rasterio.open(output_filename, 'w', **profile) as dst:
                dst.write(map_data.astype(rasterio.float32), 1)  # Cast data to float32 before writing

if __name__ == '__main__':
    create_maps()


