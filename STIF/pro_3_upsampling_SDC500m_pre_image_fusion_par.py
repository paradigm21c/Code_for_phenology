import os
import rasterio
from rasterio.enums import Resampling
from concurrent.futures import ProcessPoolExecutor

def upsample_with_rasterio(args):
    input_path, output_path, upscale_factor = args
    with rasterio.open(input_path) as dataset:
        data = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height * upscale_factor),
                int(dataset.width * upscale_factor)
            ),
            resampling=Resampling.nearest
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
     
    upscale_factor = 16
    input_directory = r'C:\Temp\MODIS\SDC500\temp' #Pre_Fusion_Crop_full'
    output_directory = r'C:\Temp\MODIS\SDC500\Pre_Fusion_Upsample_Boston'

    
    os.makedirs(output_directory, exist_ok=True)

    target_files = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith('.tif')]

    tasks = [(tgt, os.path.join(output_directory, os.path.basename(tgt).replace('stacked_', 'upsampled_')), upscale_factor) for tgt in target_files]

    with ProcessPoolExecutor(max_workers=20) as executor:
        results = list(executor.map(upsample_with_rasterio, tasks))

    print(f"SDC500m completed.")

if __name__ == '__main__':
    main()
