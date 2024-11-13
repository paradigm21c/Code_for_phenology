import os
import rasterio
from rasterio.enums import Resampling
from concurrent.futures import ProcessPoolExecutor

def upsample_with_rasterio(args):
    input_path, output_path, ref_height, ref_width = args
    with rasterio.open(input_path) as dataset:
        # Determine the scaling factors based on reference dimensions
        height_scale = ref_height / dataset.height
        width_scale = ref_width / dataset.width

        data = dataset.read(
            out_shape=(
                dataset.count,
                ref_height,
                ref_width
            ),
            resampling=Resampling.cubic_spline
        )
        transform = dataset.transform * dataset.transform.scale(
            (dataset.width / data.shape[-1]),
            (dataset.height / data.shape[-2])
        )
        with rasterio.open(output_path, 'w',
                           driver='GTiff',
                           height=data.shape[1],
                           width=data.shape[2],
                           count=dataset.count,
                           dtype=str(data.dtype),
                           crs=dataset.crs,
                           transform=transform) as dest:
            dest.write(data)

    print(f"Processed {input_path} to {output_path}")

def main():
    input_directory =  r'C:\Temp\Temperature\DayMet'
    output_directory =  r'C:\Temp\Temperature\DayMet'
    reference_directory =  r'C:\Temp\Temperature\DayMet'

    os.makedirs(output_directory, exist_ok=True)

    # Load the dimensions of the first reference file (assuming all reference files have the same dimensions)
    ref_files = [f for f in os.listdir(reference_directory) if f.endswith('tmax_20000101_30m.tif')]
    if not ref_files:
        raise ValueError("No reference files found")
    
    with rasterio.open(os.path.join(reference_directory, ref_files[0])) as ref_dataset:
        ref_height, ref_width = ref_dataset.height, ref_dataset.width

    target_files = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith('daymet.tif')]
    tasks = [(tgt, os.path.join(output_directory, os.path.basename(tgt).replace('tmax_', 'tmax2_')), ref_height, ref_width) for tgt in target_files]

    with ProcessPoolExecutor(max_workers=1) as executor:
        results = list(executor.map(upsample_with_rasterio, tasks))

    print("completed.")

if __name__ == '__main__':
    main()
