import sys
import xarray as xr
import numpy as np
import os

# Define the GFED version to use
gfed = '5'
if gfed == '4s':
    filepath = '/discover/nobackup/projects/giss_ana/users/kmezuman/gfed4s/updated'
    outputpath = '/discover/nobackup/projects/giss_ana/users/kmezuman/gfed4s/regrid'
    fname_s = 9 #number of character (starting with 0) where year starts
    fname_e = -5 #number of character (starting with -1) where year ends
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
old_nlat = 720
old_nlon = 1440
new_shape = (90, 144)  # 2x2.5 degree resolution
old_shape = (720, 1440)  # 0.25x0.25 degree resolution


# Range of years to process
years = range(1997, 2020)

# Initialize a list to store filenames and their corresponding paths
for root, dirs, files in os.walk(filepath):
    sorted_files = sorted(files) 
    nc_files_list = [os.path.join(root, file) for file in sorted_files if int(file[fname_s:fname_e]) in years]
    
# Load the latitude and longitude for the target model resolution using xarray
ds_latlon = xr.open_dataset(latpath)
if 'lat' in ds_latlon.variables and 'lon' in ds_latlon.variables:
    lat1 = ds_latlon['lat'].values
    lon1 = ds_latlon['lon'].values
else:
    raise KeyError("Variables 'lat' and 'lon' not found in the dataset.")

# Function to regrid data using xarray's interpolation method
def regrid_burned_area(burned_area_data, target_lat_points, target_lon_points):
    """
    Regrid a burned area array to a specified target resolution.
   
    Parameters:
    - burned_area_data: xarray.DataArray
        The high-resolution burned area data to be regridded.
    - target_lat_points: int
        The number of latitude points in the target resolution.
    - target_lon_points: int
        The number of longitude points in the target resolution.
       
    Returns:
    - regridded_data: xarray.DataArray
        The regridded burned area data.
    """
    # Calculate block sizes for averaging
    lat_block_size = burned_area_data.shape[0] // target_lat_points
    lon_block_size = burned_area_data.shape[1] // target_lon_points
    print(lat_block_size, lon_block_size)

    # Reshape and average over blocks
    regridded_data_values = burned_area_data.values.reshape(
            (target_lat_points, lat_block_size, target_lon_points, lon_block_size)
    ).sum(axis=(1, 3))# Average over the blocks

    # Create the target lat/lon coordinates
    lat_target = np.linspace(burned_area_data.lat.min(), burned_area_data.lat.max(), target_lat_points)
    lon_target = np.linspace(burned_area_data.lon.min(), burned_area_data.lon.max(), target_lon_points)
    print(lon_target)
    print(lon1)
    print(lat_target)
    print(lat1)

    # Create a new DataArray for the regridded data
    regridded_data = xr.DataArray(
            regridded_data_values,
            coords=[lat_target, lon_target],
            dims=["lat", "lon"]
            )
   
    return regridded_data
# Earth surface area matrix for 0.25x0.25 degree cells (in m^2)
earth_surface_area = np.full(old_shape, 510.1e12 / (old_nlon * old_nlat))  # Earth's surface area is ~510.1 million km^2

# Process each file
for filepath in nc_files_list:
    if gfed == '5':
        gfed_data = xr.open_dataset(filepath)
        # Regrid specified variables and calculate 'Nat'
        total_regridded = regrid_burned_area(gfed_data['Total'], new_nlat, new_nlon)
        if int(filepath[fname_e-4:fname_e]) > 2000:
            peat_regridded = regrid_burned_area(gfed_data['Peat'], new_nlat, new_nlon)
            crop_regridded = regrid_burned_area(gfed_data['Crop'], new_nlat, new_nlon)
            defo_regridded = regrid_burned_area(gfed_data['Defo'], new_nlat, new_nlon)
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
        gfed_data = xr.open_dataset(filepath)
        print(filepath)
        print(gfed_data.info())
        # Multiply 'small_fire_fraction' by the Earth surface area and regrid
        small_fire_fraction = gfed_data['burned_fraction']
        small_fire_area = small_fire_fraction * earth_surface_area
        small_fire_area_regridded = regrid_data_xarray(small_fire_area, lat1, lon1)
        
        # Create a new dataset for the regridded variable
        regridded_dataset = xr.Dataset({
            'small_fire_area': small_fire_area_regridded
        })
        sys.exit()
    
    # Generate unique output file name
    output_filename = os.path.join(outputpath, f"regridded_{os.path.basename(filepath)}")
    # Save the regridded dataset to a new NetCDF file
    regridded_dataset.to_netcdf(output_filename)

    print(f"Regridded data saved to {output_filename}")
