import os
import xgboost as xgb
import numpy as np
import rasterio
import re
from itertools import product
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib.pyplot as plt



def get_date_from_filename(filename):
    """
    Extracts date from the filename in format 'overlapped_YYYYMMDD_sensor.tif'
    """
    date_match = re.search(r'_(\d{8})_', filename)
    if date_match:
        return date_match.group(1)
    return None

def prepare_training_data_for_band(matched_data, band_ind, flag):
    X, y = [], []
    
    for matched_each in matched_data:
        with rasterio.open(matched_each['reference_path']) as ref_src, rasterio.open(matched_each['target_path']) as target_src:
            target_band_num_temp = []
            for bn in range(len(matched_each['band_mapping_list'])):
                # Correctly map the bands
                target_band_num_temp = matched_each['band_mapping_list'][bn]['target']
                
                if target_band_num_temp == band_ind+1:
                    band_ind = bn
                    break
                else:
                    target_band_num_temp= []
                    
            if target_band_num_temp == []:
                continue
        
            # Correctly map the bands
            ref_band_num = matched_each['band_mapping_list'][band_ind]['reference']
            target_band_num = matched_each['band_mapping_list'][band_ind]['target']
            
            # Coeff
            coeff_1 = matched_each['coeff_1'][band_ind]
            coeff_2 = matched_each['coeff_2'][band_ind]
            '''
            # Ensure that the band numbers are valid
            if ref_band_num > ref_src.count or target_band_num > target_src.count:
                continue
            '''
                
            ref_band = ref_src.read(ref_band_num)
            target_band = target_src.read(target_band_num)

            # Masks where the data is NaN for reference and target bands
            ref_nan_mask = ~np.isnan(ref_band)
            target_nan_mask = ~np.isnan(target_band)
            
            # Combine the NaN masks
            combined_nan_mask = ref_nan_mask & target_nan_mask
            
            # Combine the QA mask and the NaN masks
            valid_pixels = matched_each['qa_mask'] & combined_nan_mask
            
            X.append(target_band[valid_pixels].flatten())
            y.append(ref_band[valid_pixels].flatten())
    
    # Check if X and y are empty before attempting to concatenate
    if not X or not y:
        print("Warning: No data was appended to X or y. Skipping concatenation.")
        if flag == 'train':
            return [], [], [], []
        elif flag == 'eval':
            return [], [], [], []
            
    # Convert lists to numpy arrays and concatenate along the first dimension
    X = np.concatenate(X)
    y = np.concatenate(y)

    # Filter using valid range of Landsat ARD (7273 - 43636) and HLS (0-10000)     
    valid_mask_landsat = (X >= 7273) & (X <= 43636)
    valid_mask_hls = (y >= 0) & (y <= 10000)
    
    # Apply masks to both X and y
    X = X[valid_mask_landsat & valid_mask_hls]
    y = y[valid_mask_landsat & valid_mask_hls]
    
    # Apply scale factor
    X=X*0.0000275 -0.2
    y=y*0.0001
    
    # Transformation function (Roy et al 2016 RSE, Characterization of Landsat-7 to Landsat-8 reflective wavelength and normalized difference vegetation index continuity)
    X = X*coeff_1 + coeff_2
    
    from sklearn.model_selection import train_test_split
    X_train, X_eval, y_train, y_eval = train_test_split(X, y, test_size=0.2, random_state=42)

    if flag =='train':
        X_train, X_eval, y_train, y_eval = X_train, X_eval, y_train, y_eval
        return X_train, y_train, X_eval, y_eval
    elif flag == 'eval':
        X_train, y_train = [], []
        X_eval, y_eval = X, y
    return [], [], X_eval, y_eval

def train_xgboost(X_train, y_train, X_eval, y_eval, params):
    train_dmatrix = xgb.DMatrix(data=X_train, label=y_train)
    eval_dmatrix = xgb.DMatrix(data=X_eval, label=y_eval)
    evals = [(train_dmatrix, 'train'), (eval_dmatrix, 'eval')]
    evals_result = {}
    bst = xgb.train(params, train_dmatrix, num_boost_round=200, evals=evals, evals_result=evals_result, early_stopping_rounds=10)

    return bst

def train_band_model(matched_data, band_ind, num_cpus, num_bands, data_name):
    X_train, y_train, X_eval, y_eval = prepare_training_data_for_band(matched_data, band_ind,flag='train')
    
    params = {
        'objective': 'reg:squarederror',
        'booster': 'gbtree',
        'nthread': num_cpus // num_bands,
        'eval_metric': 'rmse',  # Change to 'rmse' for regression
    }
    print(f"Starting training for band {band_ind+1}.")
    model_band = train_xgboost(X_train.reshape(-1,1), y_train.reshape(-1,1), X_eval.reshape(-1,1), y_eval.reshape(-1,1), params)
    model_band.save_model(f"model_band_{band_ind+1}_{data_name}.json")   
    print(f"Finished training for band {band_ind+1}.")

    return model_band

def train_wrapper(args):
    # Unpack the arguments
    train_data,band_ind,num_cpus, num_bands, data_name = args
    
    return train_band_model(train_data, band_ind, num_cpus, num_bands,data_name)

def evaluate_model_on_data(model, eval_data, band_ind):
    # Prepare the data
    _, _, X_eval, y_true = prepare_training_data_for_band(eval_data, band_ind,flag='eval')
    
    # Predict using the model
    y_pred = model.predict(xgb.DMatrix(data=X_eval))
    
    # Calculate metrics
    mae = mean_absolute_error(y_true, y_pred)
    mse = mean_squared_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    
    return mae, mse, r2
    
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
    # Paths
    reference_directory = r"D:\SetoLab\HLS\Overlapped_Landsat"
    target_directory = r"D:\SetoLab\Landsat_ARD\Overlapped_SR"
    
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
    
    
    # Define clear codes for each sensor type
    clear_codes = {
        'L30': [0, 16, 32, 48, 64, 80, 96, 112, 128, 160, 176, 144],
        'S30': [0, 16, 32, 48, 64, 80, 96, 112, 128, 160, 176, 144],
        'LT04': sorted_L45789,
        'LT05': sorted_L45789,
        'LE07': sorted_L45789,
        'LC08': sorted_L45789,
        'LC09': sorted_L45789,
    }

    band_selection = {
    'L30': {'OLI':[1, 2, 3, 4, 5, 6],'others':[1, 2, 3, 4, 5, 6]},
    'S30': {'OLI':[1, 2, 3, 4, 5, 6],'others':[1, 2, 3, 4, 5, 6]},
    'LT04': [1, 2, 3, 4, 5, 6], 
    'LT05': [1, 2, 3, 4, 5, 6],
    'LE07': [1, 2, 3, 4, 5, 6],
    'LC08': [1, 2, 3, 4, 5, 6],
    'LC09': [1, 2, 3, 4, 5, 6],
    }

    coeff_1 = [0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071]
    coeff_2 = [0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]

    # List all files and get dates
    reference_files = {get_date_from_filename(f): os.path.join(reference_directory, f) for f in os.listdir(reference_directory) if f.endswith('.tif')}
    target_files = {get_date_from_filename(f): os.path.join(target_directory, f) for f in os.listdir(target_directory) if f.endswith('.tif')}

    # Match by date
    matched_files = {}
    for date, target_file in target_files.items():
        ref_file_candidates = [f for f_date, f in reference_files.items() if f_date == date]
        if not ref_file_candidates:
            continue
        
        # Prioritize the "L30" sensor if multiple candidates
        ref_file = next((f for f in ref_file_candidates if 'L30' in f), ref_file_candidates[0])
        matched_files[ref_file] = target_file

    matched_pairs = []
    for ref_file, target_file in matched_files.items():
        # Get the base filenames without extensions for output naming
        reference_dataset_path = ref_file
        target_dataset_path = target_file
    
        matched_pairs.append((reference_dataset_path, target_dataset_path))
        
    # Initialize containers for matched pair data
    matched_data = []

    for reference_dataset_path, target_dataset_path in matched_pairs:
        with rasterio.open(target_dataset_path) as tgt, rasterio.open(reference_dataset_path) as ref:
            
            sensor_name_ref = os.path.basename(reference_dataset_path).split('_')[2].split('.')[0]
            sensor_name_tgt = os.path.basename(target_dataset_path).split('_')[2].split('.')[0]
            
            # Determine which set of bands to use for reference based on target sensor
            if sensor_name_tgt in ('LC08', 'LC09'):
                ref_bands = band_selection[sensor_name_ref]['OLI']
            else:
                ref_bands = band_selection[sensor_name_ref]['others']
            
            tgt_bands = band_selection[sensor_name_tgt]
                        
            ref_clear_code = clear_codes.get(sensor_name_ref)
            tgt_clear_code = clear_codes.get(sensor_name_tgt)
            
            # Compute combined_qa_mask for the matched pair
            tgt_qa_band_data = tgt.read(tgt.count)
            ref_qa_band_data = ref.read(ref.count)
            
            tgt_clear_mask = np.isin(tgt_qa_band_data, list(tgt_clear_code))
            ref_clear_mask = np.isin(ref_qa_band_data, list(ref_clear_code))
            
            combined_qa_mask = ~tgt_clear_mask & ~ref_clear_mask
            
            band_mapping_list = [{
            'reference': ref_band, 
            'target': tgt_band
            } for ref_band, tgt_band in zip(ref_bands, tgt_bands)]  
            
            matched_data.append({
                'reference_path': reference_dataset_path,
                'target_path': target_dataset_path,
                'qa_mask': combined_qa_mask,
                'band_mapping_list': band_mapping_list,
                'coeff_1':coeff_1,
                'coeff_2':coeff_2,
            })
            
    # train_data, eval_data = train_test_split(matched_data, test_size=0.2, random_state=42)
  
    # Train a model for each band using the data from matched pairs
    models = {}
    num_cpus = 24 #cpu_count()
    num_bands = 6
    data_name = 'LandsatARD'

    with ProcessPoolExecutor(max_workers=min(num_cpus, 12)) as executor:
        for band_ind in range(num_bands):  # Corrected this loop
            
            args = (matched_data, band_ind, num_cpus, num_bands, data_name)
            
            # Submit the job to the executor and collect the future object
            executor.submit(train_band_model, *args)
   
  
if __name__ == '__main__':
    main()
    
# any code you have for loading these models is compatible with the JSON format
# loaded_model = xgb.Booster()
# loaded_model.load_model(f"model_band_{band}.json")
