import rasterio
from rasterio.transform import from_bounds
from pyproj import Proj, Transformer
import numpy as np
import concurrent.futures

# Define a function to compute lat and lon for a single row
def compute_lat_lon(row, width, utm_transform, transformer):
    lat_row = np.zeros(width, dtype=rasterio.float32)
    lon_row = np.zeros(width, dtype=rasterio.float32)
    for col in range(width):
        x, y = utm_transform * (col, row)
        lon, lat = transformer.transform(x, y)
        lat_row[col] = lat
        lon_row[col] = lon
    return row, lat_row, lon_row

# Path to your UTM18N GeoTIFF file
utm_tiff_path = r'D:\SetoLab\code\netCDF_file\Extend_PF_NYC.tif'

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
transformer = Transformer.from_proj(utm_proj, wgs84_proj)

# Parallel computation
with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = [executor.submit(compute_lat_lon, row, width, utm_transform, transformer) for row in range(height)]
    for future in concurrent.futures.as_completed(futures):
        row, lat_row, lon_row = future.result()
        lat_data[row, :] = lat_row
        lon_data[row, :] = lon_row

# Output path for the multiband GeoTIFF
output_tiff_path = r'D:\SetoLab\code\netCDF_file\Extend_PF_NYC_latlon.tif'

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
