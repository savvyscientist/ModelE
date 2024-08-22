
import sys
import xarray as xr
import numpy as np
import os
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds

# Define the GFED version to use
gfed = '5'
if gfed == '4s':
    filepath = '/discover/nobackup/projects/giss_ana/users/kmezuman/gfed4s/updated'
    outputpath = '/discover/nobackup/projects/giss_ana/users/kmezuman/gfed4s/regrid'
    fname_s = 9  # number of characters (starting with 0) where the year starts
    fname_e = -5 # number of characters (starting with -1) where the year ends
elif gfed == '5':
    filepath = '/discover/nobackup/projects/giss_ana/users/kmezuman/GFED5/BA'
    outputpath = '/discover/nobackup/projects/giss_ana/users/kmezuman/GFED5/BA_regrid2x2H'
    fname_s = 2
    fname_e = -5
else:
    raise ValueError(f"Code not compatible with GFED '{gfed}'")

# Path to the lat/lon data for the target model resolution (2x2.5 degree)
latpath = '/discover/nobackup/kmezuman/nk_CCycle_E6obioF40/ANN2005.aijnk_CCycle_E6obioF40.nc'

# Define the target resolution for regridding
new_nlon = 144
new_nlat = 90

# Range of years to process
years = range(1997, 2020)

# Initialize a list to store filenames and their corresponding paths
for root, dirs, files in os.walk(filepath):
    sorted_files = sorted(files) 
    nc_files_list = [os.path.join(root, file) for file in sorted_files if int(file[fname_s:fname_e]) in years]

# Load the latitude and longitude for the target model resolution using xarray
ds_latlon = xr.open_dataset(latpath)
if 'lat' in ds_latlon.variables and 'lon' in ds_latlon.variables:
    target_lat = ds_latlon['lat'].values
    target_lon = ds_latlon['lon'].values
else:
    raise KeyError("Variables 'lat' and 'lon' not found in the dataset.")

def read_data_with_rasterio(file_path):
    """
    Reads geospatial data using rasterio and converts it to an xarray.DataArray.

    Parameters:
    - file_path: str
        The path to the geospatial file.

    Returns:
    - data_array: xarray.DataArray
        The data as an xarray.DataArray with appropriate geospatial metadata.
    """
    with rasterio.open(file_path) as src:
        data = src.read(1)  # Read the first band
        transform = src.transform
        crs = src.crs
        bounds = src.bounds
        height, width = data.shape
        
        # Create coordinates
        lon = np.linspace(bounds.left, bounds.right, width)
        lat = np.linspace(bounds.top, bounds.bottom, height)
        
        # Create xarray DataArray
        data_array = xr.DataArray(
            data,
            coords=[lat[::-1], lon],  # Reverse latitudes to match the top-bottom convention
            dims=["lat", "lon"],
            attrs={
                'transform': transform,
                'crs': crs
            }
        )
    
    return data_array

def regrid_burned_area_rasterio(burned_area_data, target_lat, target_lon):
    """
    Regrid high-resolution burned area data to a lower-resolution grid using rasterio.

    Parameters:
    - burned_area_data: xarray.DataArray
        The high-resolution burned area data to be regridded.
    - target_lat: 
        The latitude array of the target resolution.
    - target_lon: 
        The longitude array of the target resolution.

    Returns:
    - regridded_data: xarray.DataArray
        The regridded burned area data.
    """

# Select only the first time slice if the data has a time dimension
    if 'time' in burned_area_data.dims:
        burned_area_data = burned_area_data.isel(time=0)

# Check if transform exists, otherwise create it
    if 'transform' in burned_area_data.attrs:
        src_transform = burned_area_data.attrs['transform']
    else:
        # Assuming the coordinates are 2D grids and the data is in WGS84 (EPSG:4326)
        lon_min, lon_max = burned_area_data.lon.min().item(), burned_area_data.lon.max().item()
        lat_min, lat_max = burned_area_data.lat.min().item(), burned_area_data.lat.max().item()
        src_transform = from_bounds(lon_min, lat_min, lon_max, lat_max, burned_area_data.sizes['lon'], burned_area_data.sizes['lat'])
        burned_area_data.attrs['transform'] = src_transform

# Define CRS if not present
    if 'crs' in burned_area_data.attrs:
        src_crs = burned_area_data.attrs['crs']
    else:
        src_crs = 'EPSG:4326'  # Default to WGS84
        burned_area_data.attrs['crs'] = src_crs

    
    src_height, src_width = burned_area_data.shape
    src_bounds = rasterio.transform.array_bounds(src_height, src_width, src_transform)
    
    # Define the target transform and shape
    target_transform = rasterio.transform.from_bounds(
        west=target_lon.min(),
        south=target_lat.min(),
        east=target_lon.max(),
        north=target_lat.max(),
        width=len(target_lon),
        height=len(target_lat)
    )

    target_shape = (len(target_lat), len(target_lon))

    # Initialize the destination array
    regridded_data_values = np.zeros(target_shape, dtype=burned_area_data.dtype)

    # Perform the reprojection/resampling
    reproject(
        source=burned_area_data.values,
        destination=regridded_data_values,
        src_transform=src_transform,
        src_crs=src_crs,
        src_bounds=src_bounds,
        dst_transform=target_transform,
        dst_crs=src_crs,
        resampling=Resampling.sum
    )

    # Create the xarray.DataArray for the output
    regridded_data = xr.DataArray(
        regridded_data_values,
        coords=[target_lat, target_lon],
        dims=['lat', 'lon'],
        attrs=burned_area_data.attrs
    )

    return regridded_data

def save_to_netcdf(data_array, output_path):
    """
    Save an xarray.DataArray to a NetCDF file.

    Parameters:
    - data_array: xarray.DataArray
        The data to be saved.
    - output_path: str
        The path to save the NetCDF file.
    """
    data_array.to_netcdf(output_path)
    print(f"Data successfully saved to {output_path}")

# Process each file
for filepath in nc_files_list:
    gfed_data = xr.open_dataset(filepath)

    if gfed == '5':
        # Regrid specified variables and calculate 'Nat'
        total_regridded = regrid_burned_area_rasterio(gfed_data['Total'], target_lat, target_lon)
        if int(filepath[fname_e-4:fname_e]) > 2000:
            peat_regridded = regrid_burned_area_rasterio(gfed_data['Peat'], target_lat, target_lon)
            crop_regridded = regrid_burned_area_rasterio(gfed_data['Crop'], target_lat, target_lon)
            defo_regridded = regrid_burned_area_rasterio(gfed_data['Defo'], target_lat, target_lon)
            # Calculate 'Nat' as Total - (Peat + Crop + Defo)
            nat_regridded = total_regridded - (peat_regridded + crop_regridded + defo_regridded)

            # Combine all variables into a new dataset
            regridded_dataset = xr.Dataset({
                'Total': total_regridded,
                'Peat': peat_regridded,
                'Crop': crop_regridded,
                'Defo': defo_regridded,
                'Nat': nat_regridded
            })
        else:
            regridded_dataset = xr.Dataset({
                'Total': total_regridded,
            })

    elif gfed == '4s':
        print(filepath)
        print(gfed_data.info())
        # Multiply 'small_fire_fraction' by the Earth surface area and regrid
        small_fire_fraction = gfed_data['burned_fraction']
        small_fire_area = small_fire_fraction * earth_surface_area
        small_fire_area_regridded = regrid_burned_area_rasterio(small_fire_area, target_lat, target_lon)

        # Create a new dataset for the regridded variable
        regridded_dataset = xr.Dataset({
            'small_fire_area': small_fire_area_regridded
        })
        sys.exit()

    # Generate unique output file name
    output_filename = os.path.join(outputpath, f"regridded_{os.path.basename(filepath)}")
    # Save the regridded dataset to a new NetCDF file
    save_to_netcdf(regridded_dataset, output_filename)

    print(f"Regridded data saved to {output_filename}")
