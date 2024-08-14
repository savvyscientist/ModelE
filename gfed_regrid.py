
import xarray as xr
import numpy as np
import os

# Define the GFED version to use
gfed = '4s'
if gfed == '4s':
    filepath = '/discover/nobackup/projects/giss_ana/users/kmezuman/gfed4s/updated'
elif gfed == '5':
    filepath = '/discover/nobackup/projects/giss_ana/users/kmezuman/GFED5/BA'
else:
    raise ValueError(f"Code not compatible with GFED '{gfed}'")

# Path to the lat/lon data for the target model resolution (2x2.5 degree)
latpath = '/discover/nobackup/kemzuman/ANN2005.aijnk_CCycle_E6obioF40.nc'

# Define the target resolution for regridding
new_shape = (90, 144)  # 2x2.5 degree resolution
old_shape = (720, 1440)  # 0.25x0.25 degree resolution
output_base = filepath + '/regrid'

# Initialize lists to store filenames and their corresponding paths
filename_group = []
filepath_group = []

# Range of years to process
years = range(1997, 2020)

# Generate the filenames based on the GFED version and year
for year in years:
    if gfed == '4s':
        filename = f"GFED4.1s_{year}.hdf5"
    elif gfed == '5':
        for x in range(12):
            if x < 9:
                filename = f"BA{year}0{x + 1}.nc"
            else:           
                filename = f"BA{year}{x + 1}.nc"
    else:
        raise ValueError(f"Code not compatible with GFED '{gfed}'")
    filename_group.append(filename)

# Generate the full file paths for the data
for m in filename_group:
    filepath_place = os.path.join(filepath, m)
    filepath_group.append(filepath_place)

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

# Earth surface area matrix for 0.25x0.25 degree cells (in km^2)
earth_surface_area = np.full(old_shape, 510.1e6 / (720 * 1440))  # Earth's surface area is ~510.1 million km^2

# Process each file
for filepath in filepath_group:
    # Open the GFED dataset
    gfed_data = xr.open_dataset(filepath)
    
    if gfed == '5':
        # Regrid specified variables and calculate 'Nat'
        total_regridded = regrid_data_xarray(gfed_data['Total'], lat1, lon1)
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
        
    elif gfed == '4s':
        # Multiply 'small_fire_fraction' by the Earth surface area and regrid
        small_fire_fraction = gfed_data['small_fire_fraction']
        small_fire_area = small_fire_fraction * earth_surface_area
        small_fire_area_regridded = regrid_data_xarray(small_fire_area, lat1, lon1)
        
        # Create a new dataset for the regridded variable
        regridded_dataset = xr.Dataset({
            'small_fire_area': small_fire_area_regridded
        })
    
    # Save the regridded dataset to a new NetCDF file
    output_filepath = os.path.join(output_base, os.path.basename(filepath).replace('.hdf5', '_regridded.nc').replace('.nc', '_regridded.nc'))
    regridded_dataset.to_netcdf(output_filepath)

    print(f"Regridded data saved to {output_filepath}")
