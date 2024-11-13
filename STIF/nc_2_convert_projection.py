import os
from osgeo import gdal, osr
import glob

gdal.UseExceptions()

def reproject_tiff(input_tiff, output_tiff, src_crs, tgt_crs):
    gdal.UseExceptions()
    try:
        # Open the source file
        src_ds = gdal.Open(input_tiff)

        # Reproject
        reprojected_ds = gdal.Warp(output_tiff, src_ds, srcSRS=src_crs, dstSRS=tgt_crs)

        # Clean up
        del reprojected_ds, src_ds
    except RuntimeError as e:
        print(f"An error occurred: {e}")


# Define source CRS (LCC)
source_crs = osr.SpatialReference()
source_crs.ImportFromProj4(
    '+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')

# Define target CRS (UTM Zone 18N)
target_crs = osr.SpatialReference()
target_crs.ImportFromEPSG(32618)

input_folder = 'C:\\Temp\\Temperature\\DayMet2\\'
output_folder = 'C:\\Temp\\Temperature\\DayMet2\\DayMet2_UTM18N\\'

for input_tiff in glob.glob(input_folder + '*.tif'):
    output_tiff = output_folder + os.path.basename(input_tiff)
    reproject_tiff(input_tiff, output_tiff, source_crs, target_crs)