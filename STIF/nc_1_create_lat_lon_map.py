import rasterio
from rasterio.transform import from_bounds
from pyproj import Proj, transform
import numpy as np

# Path to your UTM18N GeoTIFF file
utm_tiff_path = r'D:\SetoLab\code\netCDF_file\Extend_MODIS_NYC.tif'

# Open the GeoTIFF file
with rasterio.open(utm_tiff_path) as utm_tiff:
    utm_transform = utm_tiff.transform
    utm_crs = utm_tiff.crs
    width = utm_tiff.width
    height = utm_tiff.height

# Initialize latitude and longitude arrays
lat_data = np.zeros((height, width), dtype=rasterio.float32)
lon_data = np.zeros((height, width), dtype=rasterio.float32)

# Set up projections
utm_proj = Proj(utm_crs)  # UTM Zone 18N projection
wgs84_proj = Proj(proj='latlong', datum='WGS84')

# Compute latitude and longitude for each pixel
for row in range(height):
    for col in range(width):
        x, y = utm_transform * (col, row)
        lon, lat = transform(utm_proj, wgs84_proj, x, y)
        lat_data[row, col] = lat
        lon_data[row, col] = lon

# Output path for the multiband GeoTIFF
output_tiff_path = r'D:\SetoLab\code\netCDF_file\Extend_MODIS_NYC_latlon.tif'

# Create a new multiband GeoTIFF
with rasterio.open(
    output_tiff_path, 
    'w', 
    driver='GTiff',
    height=height, 
    width=width,
    count=2,  # Two bands: latitude and longitude
    dtype=rasterio.float32,
    crs=utm_crs,  # Keep the same CRS
    transform=utm_transform
) as dst:
    dst.write(lat_data, 1)  # Write latitude data to band 1
    dst.write(lon_data, 2)  # Write longitude data to band 2
