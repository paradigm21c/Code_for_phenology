import os
import numpy as np
import rasterio
from concurrent.futures import ProcessPoolExecutor
from itertools import product

def check_clearness(dataset_path, clear_codes):
    with rasterio.open(dataset_path) as dataset:
        # Read the last band (QA band)
        qa_band = dataset.read(dataset.count)
        
        # Check for nan values
        nan_mask = np.isnan(qa_band)
        
        # Check for non-clear pixels using clear codes
        clear_mask = np.isin(qa_band, clear_codes)
        
        fill_mask = np.isin(qa_band, 1)
        
        # Combine the masks: True for contaminated pixels, False for clear pixels
        contaminated_mask = ~clear_mask | nan_mask | fill_mask
        
        # Calculate percentage of contaminated pixels
        contaminated_percentage = np.sum(contaminated_mask) / (dataset.width * dataset.height) * 100
        
        return contaminated_percentage < 20 , contaminated_mask

def process_single_file(filepath, output_directory, clear_codes_map):
    sensor_type = os.path.basename(filepath).split('_')[2].split('.')[0]
    clear_codes = clear_codes_map.get(sensor_type)
    
    is_clear, mask = check_clearness(filepath, clear_codes)
    
    if is_clear:
        # Save with a different name if it's a clear image
        output_path = os.path.join(output_directory, os.path.basename(filepath).replace('cropped','clear'))
        with rasterio.open(filepath) as src:
            # Read all bands except the last one (QA band)
            bands_data = src.read(indexes=list(range(1, src.count)))
            
            # Ensure the data is in floating point format to handle np.nan
            bands_data = bands_data.astype('float32')
            
            # Expand the mask to match the number of bands
            # expanded_mask = np.repeat(mask[np.newaxis, :, :], bands_data.shape[0], axis=0)
            
            # Apply the mask
            # bands_data[expanded_mask] = np.nan
            
            # Adjust metadata for the number of bands and data type
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

def main():

    # Specify known bit conditions 
    known_bits = {
        15: ['0'],  # not fill (LSB)
        12: ['0'],  # Cloud
        11: ['0'],  # Cloud shadow
        (7, 6): ['00', '01', '10'],  # none and small cirrus
        (5, 4): ['00', '01', '10']   # conditions for bits 3 and 4
    }

    # For the unspecified bits, allow any combination
    unspecified_bits_indices = [i for i in range(16) if i not in known_bits and i not in [j for k in known_bits.keys() if isinstance(k, tuple) for j in range(k[0], k[1] + 1)]]
    for idx in unspecified_bits_indices:
        known_bits[idx] = ['0', '1']

    # Expand the ranges into individual bit positions
    expanded_bits = expand_ranges(known_bits)

    # Generate all possible combinations considering the order of bits
    sorted_keys = sorted(expanded_bits.keys())
    combinations = list(product(*[expanded_bits[key] for key in sorted_keys]))

    # Convert bit strings to integers
    int_values = [int(''.join(combination), 2) for combination in combinations]
    sorted_L45789 = sorted(list(set(int_values)))
    

    ## Define clear codes for each sensor type   
    clear_codes_map = {
        'LT04': sorted_L45789,
        'LT05': sorted_L45789,
        'LE07': sorted_L45789,
        'LC08': sorted_L45789,
        'LC09': sorted_L45789,
    }

    input_directory = r"D:\SetoLab\Landsat_ARD\Pre_Fusion_Crop_full_ST" 
    output_directory = r"D:\SetoLab\Landsat_ARD\Pre_Fusion_Clear_full_ST"  

    # Ensure output directories exist
    os.makedirs(output_directory, exist_ok=True)

    # Get list of all files to process
    filepaths = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith('.tif')]

    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=20) as executor:
        executor.map(process_single_file, filepaths, [output_directory]*len(filepaths), [clear_codes_map]*len(filepaths))

if __name__ == "__main__":
    main()