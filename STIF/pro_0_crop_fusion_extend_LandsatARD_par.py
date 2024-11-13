import os
import re
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import Window
from datetime import datetime, timedelta
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor


def get_overlapping_window(src1, src2):
    """
    Determines the overlapping window of two raster datasets.
    Returns the window for both rasters as a tuple.
    """
    # Get the bounds of both datasets
    left1, bottom1, right1, top1 = src1.bounds
    left2, bottom2, right2, top2 = src2.bounds
    
    # Determine the overlapping bounds
    left, bottom, right, top = max(left1, left2), max(bottom1, bottom2), min(right1, right2), min(top1, top2)
    
    # If there's no overlap, return None
    if left >= right or bottom >= top:
        return None, None

    # Get the window of the overlap for each dataset
    window1 = src1.window(left, bottom, right, top)
    window2 = src2.window(left, bottom, right, top)
    
    return window1, window2

def extract_overlap_from_datasets(reference_path, target_path, target_out_path):
    """Extracts overlapping regions from two raster datasets."""
    with rasterio.open(reference_path) as ref_src, rasterio.open(target_path) as tgt_src:
        
        window1, window2 = get_overlapping_window(ref_src, tgt_src)

        # If there's no overlap, skip further processing
        if window1 is None or window2 is None:
            print(f"No overlap found between {reference_path} and {target_path}.")
            return

        # Read the overlapping region data directly without resampling
        tgt_data = tgt_src.read(window=window2)

        # Update metadata        
        meta2 = tgt_src.meta.copy()
        meta2.update({
            'width': window2.width,
            'height': window2.height,
            'transform': rasterio.windows.transform(window2, tgt_src.transform)
        })

        # Save the overlapping regions
        with rasterio.open(target_out_path, 'w', **meta2) as tgt_out:
            tgt_out.write(tgt_data)
            
def extract_date_from_filename(filename):
    """Extract date from filename using regex."""
    match = re.search(r"(\d{8})", filename)
    if match:
        return datetime.strptime(match.group(1), '%Y%m%d')
    return None

def process_target_file(reference_file, tgt_file, output_directory_target):
    tgt_out_name = os.path.basename(tgt_file).replace('reprojected', 'cropped')
    tgt_out_path = os.path.join(output_directory_target, tgt_out_name)
    # Extract overlapping areas
    extract_overlap_from_datasets(reference_file, tgt_file, tgt_out_path)
    print(f"Processed prefusion for {tgt_file}.")
    
def main():
    # Paths
    reference_directory = r"D:\SetoLab\code\Processing"
    target_directory = r"C:\Temp\Landsat_ARD\Reprojected_Extracted_ST"
    output_directory_target = r"D:\SetoLab\Landsat_ARD\Pre_Fusion_Crop_full_ST"

    # Ensure output directories exist
    os.makedirs(output_directory_target, exist_ok=True)

    # List of files
    reference_files = [os.path.join(reference_directory, f) for f in os.listdir(reference_directory) if f.endswith('Fusion_NYC_full.tif')]
    target_files = [os.path.join(target_directory, f) for f in os.listdir(target_directory) if f.endswith('.tif')]

    # Assuming there's only one reference file
    if len(reference_files) != 1:
        raise ValueError("There should be only one reference file.")

    reference_file = reference_files[0] 
    num_cpus = cpu_count()
    
    # Parallelize the processing of target files using a ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        list(executor.map(process_target_file, 
                         [reference_file] * len(target_files), 
                         target_files, 
                         [output_directory_target] * len(target_files)))

if __name__ == '__main__':
    main()
