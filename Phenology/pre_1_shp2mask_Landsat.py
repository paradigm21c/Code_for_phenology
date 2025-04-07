import os
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from concurrent.futures import ProcessPoolExecutor

def create_mask(shp_file, geotiff_path, maskfile_path, shapefile_path):
    shapes = gpd.read_file(os.path.join(shapefile_path, shp_file))

    with rasterio.open(geotiff_path) as src:
        out_image, out_transform = mask(src, shapes.geometry, all_touched=False, crop=False, indexes=1)
        out_meta = src.meta
        out_image[out_image !=src.meta['nodata']]=1
        out_image[out_image ==src.meta['nodata']]=0

    out_meta.update({"driver": "GTiff",
                     "count": 1,
                     "height": src.height,
                     "width": src.width,
                     "nodata": -9999,
                     "transform": src.transform})

    output_name = shp_file.replace('.shp', '_Landsat.tif')
    with rasterio.open(os.path.join(maskfile_path, output_name), 'w', **out_meta) as dest:
        dest.write(out_image, 1)

def main():
    maskfile_path = 'D:/SetoLab/Phenology/mask'
    shapefile_path = 'D:/SetoLab/NAC_ECM/SHP'
    geotiff_path = 'D:/SetoLab/One30m/Pre_Fusion_500m/clear_20230902_LC09.tif' # Landsat

    shp_files = [file for file in os.listdir(shapefile_path) if file.endswith('.shp')]

    # Using ProcessPoolExecutor to parallelize the task
    with ProcessPoolExecutor(max_workers=5) as executor:
        for shp_file in shp_files:
            executor.submit(create_mask, shp_file, geotiff_path, maskfile_path, shapefile_path)

if __name__ == "__main__":
    main()
