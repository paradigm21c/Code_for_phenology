import os
import rasterio
import numpy as np
from sklearn.cluster import KMeans
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor

def classify_image(input_path, output_dir, n_clusters):
    # Read the raster data
    with rasterio.open(input_path) as src:
        data = src.read()
        meta = src.meta

    # Get the shape of the data
    bands, height, width = data.shape

    # Reshape the data to a 2D array, where each row is a pixel and each column is a band
    reshaped_data = data.reshape(bands, -1).T
    
    # Set values outside the [0, 10000] range to NaN
    mask = np.logical_or(reshaped_data < 0, reshaped_data > 10000)
    reshaped_data[mask] = 0

    # Apply KMeans clustering
    kmeans = KMeans(n_clusters=n_clusters, n_init='auto')
    kmeans_labels = kmeans.fit_predict(reshaped_data)

    # Reshape the clustered data back to the original 2D shape
    classified_data = kmeans_labels.reshape(height, width)

    # Update metadata to save as single band
    meta.update(count=1, dtype='int16', nodata=-9999)
    
    out_name = os.path.basename(input_path).replace('clear', 'classified')       
    out_path = os.path.join(output_dir, out_name)

    # Write the classified data to a new raster file
    with rasterio.open(out_path, 'w', **meta) as dst:
        dst.write(classified_data.astype(rasterio.int16), 1)
        
    return out_path
    
if __name__ == '__main__':
    input_directory = r"D:/SetoLab/One30m/Pre_Fusion_500m_full"
    output_directory = r"D:/SetoLab/One30m/Pre_Fusion_Class_full"

    # Ensure output directories exist
    os.makedirs(output_directory, exist_ok=True)
        
    # Get list of all files to process
    files = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith('.tif')]     
    number_of_classes = 8

    # Parallelize using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=20) as executor:
        results = list(executor.map(classify_image, files, [output_directory]*len(files), [number_of_classes]*len(files)))

    for r in results:
        print(f"Classified {r}")