import os
import xgboost as xgb
import rasterio
from datetime import datetime
from itertools import repeat
import numpy as np
from concurrent.futures import ProcessPoolExecutor

def apply_model_to_band(band_data, model_path):
    """
    Apply the XGBoost model to the given band data.
    """
    # Reshape band data for prediction
    original_shape = band_data.shape
    flattened_data = band_data.reshape(-1, 1)  # Convert 2D array into long 2D array with one feature
    
    dmatrix = xgb.DMatrix(flattened_data)
    model = xgb.Booster()
    model.load_model(model_path)
    predictions = model.predict(dmatrix)
    
    # Reshape the predictions back to the original band shape
    harmonized_band_data = predictions.reshape(original_shape)
    
    return harmonized_band_data


def harmonize_dataset(dataset_path, output_directory, coeffs):
    with rasterio.open(dataset_path) as src:
        meta = src.meta
        harmonized_data = []
        for band_num in range(1, meta['count']+1):
            band_data = src.read(band_num)
            band_data = band_data*0.0000275 - 0.2
            
            # Apply harmonization to all but the last band
            if band_num < meta['count']:
                coeff_1 = coeffs['coeff_1'][band_num - 1]
                coeff_2 = coeffs['coeff_2'][band_num - 1]
                band_data = band_data * coeff_1 + coeff_2

                model_path = f"model_band_{band_num}_LandsatARD.json"
                harmonized_band_data = apply_model_to_band(band_data, model_path)
                harmonized_band_data = harmonized_band_data * 10000
                harmonized_band_data = harmonized_band_data.astype(np.int16)  # Correct data type conversion
                harmonized_data.append(harmonized_band_data)
            else:
                # For the last band, just append it without harmonization
                last_band_data = src.read(meta['count'])
                last_band_data = last_band_data.astype(np.int16)
                harmonized_data.append(last_band_data)
                
        # Construct the output filename
        sensor_name = os.path.basename(dataset_path).split('_')[0]
        output_filename = os.path.basename(dataset_path).replace('cropped', 'harmonized')
        output_path = os.path.join(output_directory, output_filename)

        # Update meta to reflect the number of bands and dtype
        meta.update(count=len(harmonized_data), dtype='int16')
        
        # Save harmonized data
        with rasterio.open(output_path, 'w', **meta) as dst:
            for band_num, band in enumerate(harmonized_data, 1):
                dst.write(band, band_num)

    return f"Processed {dataset_path}"

if __name__ == '__main__':
    INPUT_DIR = r"D:\SetoLab\Landsat_ARD\Pre_Fusion_Crop_full"
    OUTPUT_DIR = r"D:\SetoLab\Landsat_ARD\Pre_Fusion_Harmonized_full"
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    coeffs = {
    'coeff_1': [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071],
    'coeff_2': [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]
    }
    
    dataset_paths = [os.path.join(INPUT_DIR, f) for f in os.listdir(INPUT_DIR) if f.endswith('.tif')]
    
    # Parallel processing using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=12) as executor:
        results = list(executor.map(harmonize_dataset, dataset_paths, repeat(OUTPUT_DIR), repeat(coeffs)))
    
    for result in results:
        print(result)
