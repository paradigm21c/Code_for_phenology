import rasterio
import numpy as np
import xarray as xr
from rasterio.transform import from_origin
from scipy.spatial import cKDTree
import datetime
import concurrent.futures
import os

def daymet_index_to_date(year, day_index):
    # Daymet treats leap years by including Feb 29 and excluding Dec 31.
    # Thus, for leap years, day 60 (Feb 29) exists, but day 366 (Dec 31) does not.
    # Simply return the date, assuming a 365-day year
    return datetime.date(year, 1, 1) + datetime.timedelta(days=day_index - 1)
    
# Function to find nearest index
def find_nearest_index(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def compute_nearest_tmax_row(y, lon_band_row, lat_band_row, tree, tmax_data, width):
    nearest_tmax_row = np.zeros(width, dtype=np.float32)
    for x in range(width):
        lat = lat_band_row[x]
        lon = lon_band_row[x]
        _, idx = tree.query([lat, lon])
        nearest_tmax_row[x] = tmax_data.ravel()[idx]
    return y, nearest_tmax_row
    



product = 'MODIS'
    
# Path to your latitude and longitude GeoTIFF file
lat_lon_tiff_path = f'D:/SetoLab/code/netCDF_file/Extend_{product}_NYC_latlon.tif'

# Read the latitude and longitude bands
with rasterio.open(lat_lon_tiff_path) as lat_lon_tiff:
    lat_band = lat_lon_tiff.read(1)
    lon_band = lat_lon_tiff.read(2)
    tiff_transform = lat_lon_tiff.transform
    tiff_crs = lat_lon_tiff.crs


for vari in ['vp']: #,'srad']: # 'prcp',

    for yr in range(1999,2000):
        # Path to your NetCDF file
        nc_file = f'C:\Temp\Temperature\DayMet\raw\daymet_v4_daily_na_{vari}_{yr}.nc' #f'H:\DayMet_preci\daymet_v4_daily_na_{vari}_{yr}.nc'
        # ds = xr.open_dataset(nc_file)
        ds = xr.open_dataset(nc_file, chunks={"time": 1})

        # Extract the 'tmax' variable
        tmax = ds[f'{vari}']
        lats = ds['lat'].values
        lons = ds['lon'].values

        # Create a KDTree for quick nearest-neighbor lookup
        lat_lon_pairs = np.column_stack([lats.ravel(), lons.ravel()])
        tree = cKDTree(lat_lon_pairs)

        # Loop through each time step and save as a GeoTIFF
        for i, timestep in enumerate(tmax.time):
            # Extract data for the day
            tmax_data = tmax.isel(time=i).values
            
            year = timestep.dt.year.item()  # Directly get the year from the time_step
            date = daymet_index_to_date(year, i + 1)
            date_str = date.strftime('%Y%m%d')

            # Initialize an array to hold nearest 'tmax' values
            nearest_tmax_data = np.full(lat_band.shape, np.nan, dtype=np.float32)

            # Parallel computation
            with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
                futures = []
                for y in range(lat_band.shape[0]):
                    future = executor.submit(compute_nearest_tmax_row, y, lon_band[y], lat_band[y], tree, tmax_data, lat_band.shape[1])
                    futures.append(future)

                for future in concurrent.futures.as_completed(futures):
                    y, nearest_tmax_row = future.result()
                    nearest_tmax_data[y, :] = nearest_tmax_row
            
            output_directory_target = f'C:\\Temp\\Temperature\\DayMet\\{product}\\{vari}\\{yr}\\'
            os.makedirs(output_directory_target, exist_ok=True)
            # Define output GeoTIFF path
            output_tiff_path = os.path.join(output_directory_target,f'{vari}_{date_str}_daymet.tif')

            # Create and save the GeoTIFF
            with rasterio.open(
                output_tiff_path, 
                'w',
                driver='GTiff',
                height=nearest_tmax_data.shape[0],
                width=nearest_tmax_data.shape[1],
                count=1,
                dtype=nearest_tmax_data.dtype,
                crs=tiff_crs,
                transform=tiff_transform
            ) as dst:
                dst.write(nearest_tmax_data, 1)
