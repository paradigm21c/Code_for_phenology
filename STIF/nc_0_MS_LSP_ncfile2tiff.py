import xarray as xr
import rasterio
from rasterio.transform import Affine
import numpy as np

path = 'C:/Temp/HLS/MuSLI_LSP'


# Function to convert DataArray to GeoTIFF with UTM18N projection
def convert_to_geotiff(path, data_array, variable_name, yr, tile, tm):
    # Convert the xarray DataArray to a numpy array
    np_array = data_array.values

    # Extract the GeoTransform parameters from the attribute
    geo_transform = tm.attrs['GeoTransform'].split()
    geo_transform = [float(x) for x in geo_transform]
    transform = Affine.from_gdal(*geo_transform)

    # Define the output GeoTIFF file path
    tiff_file = f'{path}/{tile}_{variable_name}_{yr}.tif'

    # Define metadata
    metadata = {
        'driver': 'GTiff',
        'height': np_array.shape[0],  # Assuming the shape order is (y, x)
        'width': np_array.shape[1],  # Assuming the shape order is ( y, x)
        'count': 1,
        'dtype': np_array.dtype,
        'crs': 'EPSG:32618',  # UTM Zone 18N
        'transform': transform,
    }

    # Write the GeoTIFF file
    with rasterio.open(tiff_file, 'w', **metadata) as dst:
        dst.write(np_array[:, :], 1)  # Write the first time step

    print(f"GeoTIFF file for variable '{variable_name}' saved as {tiff_file}")

# Iterate over each year and convert variables to GeoTIFF
tiles = ['18TWK', '18TWL', '18TXL']

for tile in tiles:
    for yr in range(2016, 2017):
        nc_file = f'{path}/MSLSP_{tile}_{yr}.nc'  # Update the file path to include the year
        ds = xr.open_dataset(nc_file, engine='netcdf4')
        variables = list(ds.data_vars)
        tm = ds['transverse_mercator']  # Extract the transverse_mercator variable for GeoTransform
        variables = variables[1:]  # Skip the transverse_mercator variable
        print(f"Variables in the NetCDF file for year {yr}: {variables}")

        for variable in variables:
            data_array = ds[variable]
            convert_to_geotiff(path, data_array, variable, yr, tile, tm)
