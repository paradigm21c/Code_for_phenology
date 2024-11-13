import netCDF4 as nc
import xarray as xr
import numpy as np
import rasterio
from rasterio.transform import from_origin
import datetime

# # Path to your NetCDF file
# file_path = 'C:\Temp\Temperature\DayMet2\daymet_v4_daily_na_tmax_2016.nc'
#
# # Using netCDF4
# ds_nc = nc.Dataset(file_path)
# print("Using netCDF4:")
# print("Variables:", ds_nc.variables.keys())
# print("Dimensions:", ds_nc.dimensions.keys())
# for var in ds_nc.variables.values():
#     print(var)  # Print details of each variable
#
# # Using xarray
# ds_xr = xr.open_dataset(file_path)
# print("\nUsing xarray:")
# print(ds_xr.info())
#

def daymet_index_to_date(year, day_index):
    # Daymet treats leap years by including Feb 29 and excluding Dec 31.
    # Thus, for leap years, day 60 (Feb 29) exists, but day 366 (Dec 31) does not.
    # Simply return the date, assuming a 365-day year
    return datetime.date(year, 1, 1) + datetime.timedelta(days=day_index - 1)

def save_to_geotiff(data_slice, output_filename, transform, crs):
    with rasterio.open(
        output_filename,
        'w',
        driver='GTiff',
        height=data_slice.shape[0],
        width=data_slice.shape[1],
        count=1,
        dtype=data_slice.dtype,
        crs=crs,
        transform=transform,
    ) as dst:
        dst.write(data_slice, 1)


# Open the NetCDF file
nc_file = 'C:\Temp\Temperature\DayMet2\daymet_v4_daily_na_tmax_2016.nc'
ds = xr.open_dataset(nc_file)

# Extract the 'tmax' variable
tmax = ds['tmax']

# Calculate resolutions for x and y
x_resolution = ds['x'][1] - ds['x'][0]
y_resolution = ds['y'][1] - ds['y'][0]

# It's also a good practice to check if the resolution is consistent throughout the dataset
assert np.allclose(np.diff(ds['x']), x_resolution, atol=1e-6)
assert np.allclose(np.diff(ds['y']), y_resolution, atol=1e-6)

# Define CRS and transform for Lambert Conformal Conic
crs = "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

# Define the transform
transform = from_origin(min(ds['x'].values), max(ds['y'].values), x_resolution, y_resolution)

start_year = 2016  # Adjust this to your dataset's start year

for i, time_step in enumerate(tmax.time):
    year = time_step.dt.year.item()  # Directly get the year from the time_step
    date = daymet_index_to_date(year, i + 1)
    date_str = date.strftime('%Y%m%d')

    data_slice = tmax.isel(time=i).values
    output_filename = f'C:\\Temp\\Temperature\\DayMet2\\tmax_{date_str}_daymet.tif'
    save_to_geotiff(data_slice, output_filename, transform, crs)
