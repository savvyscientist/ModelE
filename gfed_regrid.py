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
def regrid_data_xarray(data, new_lat, new_lon):
    regridded = data.interp(lat=new_lat, lon=new_lon, method='linear')
    return regridded

# Earth surface area matrix for 0.25x0.25 degree cells (in m^2)
earth_surface_area = np.full(old_shape, 510.1e12 / (old_nlon * old_nlat))  # Earth's surface area is ~510.1 million km^2

# Process each file
for filepath in nc_files_list:
    # Open the GFED dataset
    gfed_data = xr.open_dataset(filepath)
    print(filepath)
    #print(gfed_data.info())
    
    if gfed == '5':
        # Regrid specified variables and calculate 'Nat'
        total_regridded = regrid_data_xarray(gfed_data['Total'], lat1, lon1)
        if int(filepath[fname_e-4:fname_e]) > 2000:
            peat_regridded = regrid_data_xarray(gfed_data['Peat'], lat1, lon1)
            crop_regridded = regrid_data_xarray(gfed_data['Crop'], lat1, lon1)
            defo_regridded = regrid_data_xarray(gfed_data['Defo'], lat1, lon1)
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
        # Multiply 'small_fire_fraction' by the Earth surface area and regrid
        small_fire_fraction = gfed_data['burned_fraction']
        small_fire_area = small_fire_fraction * earth_surface_area
        small_fire_area_regridded = regrid_data_xarray(small_fire_area, lat1, lon1)
        
        # Create a new dataset for the regridded variable
        regridded_dataset = xr.Dataset({
            'small_fire_area': small_fire_area_regridded
        })
    
    # Save the regridded dataset to a new NetCDF file
    regridded_dataset.to_netcdf(outputpath)

    print(f"Regridded data saved to {outputpath}")
