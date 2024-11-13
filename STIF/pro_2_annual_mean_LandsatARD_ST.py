import os
import numpy as np
import rasterio
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from collections import defaultdict

def check_clearness(dataset_path, clear_codes):
    with rasterio.open(dataset_path) as dataset:
        qa_band = dataset.read(dataset.count)
        nan_mask = np.isnan(qa_band)
        clear_mask = np.isin(qa_band, clear_codes)
        fill_mask = np.isin(qa_band, 1)
        contaminated_mask = ~clear_mask | nan_mask | fill_mask
        contaminated_percentage = np.sum(contaminated_mask) / (dataset.width * dataset.height) * 100
        return contaminated_percentage < 20, contaminated_mask

def process_single_file(filepath, output_directory, clear_codes_map):
    sensor_type = os.path.basename(filepath).split('_')[2].split('.')[0]
    clear_codes = clear_codes_map.get(sensor_type)
    is_clear, mask = check_clearness(filepath, clear_codes)
    if is_clear:
        output_path = os.path.join(output_directory, os.path.basename(filepath).replace('cropped', 'clear'))
        with rasterio.open(filepath) as src:
            bands_data = src.read(indexes=list(range(1, src.count)))
            bands_data = bands_data.astype('float32')
            bands_data = bands_data*0.00341802 + 149.0
            
            for band_index in range(bands_data.shape[0]):
                bands_data[band_index][mask] = np.nan
            meta = src.meta
            meta['count'] = src.count - 1
            meta['dtype'] = 'float32'
            with rasterio.open(output_path, 'w', **meta) as dst:
                dst.write(bands_data)

def expand_ranges(known_bits):
    expanded_bits = {}
    for key, values in known_bits.items():
        if isinstance(key, tuple):
            for i in range(key[0], key[1] + 1):
                expanded_bits[i] = values
        else:
            expanded_bits[key] = values
    return expanded_bits

def calculate_annual_mean(files, annual_mean_output_path):
    with rasterio.open(files[0]) as src:
        meta = src.meta.copy()
        meta.update(dtype='float32')
        annual_mean = np.zeros((meta['count'], meta['height'], meta['width']), dtype='float32')
        valid_pixel_count = np.zeros((meta['count'], meta['height'], meta['width']), dtype='float32')
    
    for file in files:
        with rasterio.open(file) as src:
            data = src.read()
            valid_mask = ~np.isnan(data)
            annual_mean += np.where(valid_mask, data, 0)
            valid_pixel_count += valid_mask.astype(float)
    
    nonzero_mask = valid_pixel_count > 0
    annual_mean[nonzero_mask] /= valid_pixel_count[nonzero_mask]
    annual_mean[~nonzero_mask] = np.nan

    with rasterio.open(annual_mean_output_path, 'w', **meta) as dst:
        dst.write(annual_mean)

def main():
    known_bits = {
        15: ['0'], 12: ['0'], 11: ['0'], (7, 6): ['00', '01', '10'], (5, 4): ['00', '01', '10']
    }
    unspecified_bits_indices = [i for i in range(16) if i not in known_bits and i not in [j for k in known_bits.keys() if isinstance(k, tuple) for j in range(k[0], k[1] + 1)]]
    for idx in unspecified_bits_indices:
        known_bits[idx] = ['0', '1']
    expanded_bits = expand_ranges(known_bits)
    sorted_keys = sorted(expanded_bits.keys())
    combinations = list(product(*[expanded_bits[key] for key in sorted_keys]))
    int_values = [int(''.join(combination), 2) for combination in combinations]
    sorted_L45789 = sorted(list(set(int_values)))

    clear_codes_map = {
        'LT04': sorted_L45789, 'LT05': sorted_L45789, 'LE07': sorted_L45789, 'LC08': sorted_L45789, 'LC09': sorted_L45789,
    }

    input_directory = r"D:\SetoLab\Landsat_ARD\Pre_Fusion_Crop_ST" 
    output_directory = r"D:\SetoLab\Landsat_ARD\Pre_Fusion_Clear_ST"
    annual_mean_output_directory = r"D:\SetoLab\Landsat_ARD\Annual_Mean"

    os.makedirs(output_directory, exist_ok=True)
    os.makedirs(annual_mean_output_directory, exist_ok=True)

    filepaths = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith('.tif')]

    with ProcessPoolExecutor(max_workers=20) as executor:
        executor.map(process_single_file, filepaths, [output_directory]*len(filepaths), [clear_codes_map]*len(filepaths))

    clear_files = [os.path.join(output_directory, f) for f in os.listdir(output_directory) if f.startswith('clear')]

    files_by_year_sensor = defaultdict(list)
    for file in clear_files:
        filename = os.path.basename(file)
        year = filename.split('_')[1][:4]
        sensor = filename.split('_')[2].split('.')[0]
        files_by_year_sensor[(year, sensor)].append(file)

    for (year, sensor), files in files_by_year_sensor.items():
        annual_mean_output_path = os.path.join(annual_mean_output_directory, f'mean_{year}_{sensor}.tif')
        calculate_annual_mean(files, annual_mean_output_path)
        

if __name__ == "__main__":
    main()
