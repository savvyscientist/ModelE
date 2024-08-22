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
#old_nlat = 720
#old_nlon = 1440
#old_shape = (720, 1440)  # 0.25x0.25 degree resolution


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

# Function to regrid data using xarray's interpolation method
def regrid_burned_area(burned_area_data, target_lat, target_lon):
    """
    Manually perform conservative regridding from a high-resolution grid to a lower-resolution grid.
    
    Parameters:
    - burned_area_data: xarray.DataArray
        The high-resolution burned area data to be regridded.
    - target_lat: 
        The latitude array of the target resolution.
    - target_lon_points: int
        The longitude array of the target resolution.
        
    Returns:
    - regridded_data: xarray.DataArray
        The regridded burned area data.
    """
    # Read the original grid
    source_lat = burned_area_data.lat
    source_lon = burned_area_data.lon
    
    # Initialize the new regridded array
    regridded_data_values = np.zeros((target_lat.size, target_lon.size))
    
    # Calculate the grid cell size for both grids
    lat_old_cell_size = 180./source_lat.size
    lon_old_cell_size = 360./source_lon.size
    lat_new_cell_size = 180./target_lat.size
    lon_new_cell_size = 360./target_lon.size 
    # Loop over each cell in the target grid

    for i in range(target_lat):
        for j in range(target_lon):
            # Find the bounds of the new grid cell
            lat_min = lat_new[i] - lat_new_cell_size / 2
            lat_max = lat_new[i] + lat_new_cell_size / 2
            lon_min = lon_new[j] - lon_new_cell_size / 2
            lon_max = lon_new[j] + lon_new_cell_size / 2

            # Find which cells in the original grid contribute to this new cell
            lat_indices = np.where((lat_old >= lat_min) & (lat_old < lat_max))[0]
            lon_indices = np.where((lon_old >= lon_min) & (lon_old < lon_max))[0]

            # Accumulate the values from the original grid into the new grid cell
            for lat_index in lat_indices:
                for lon_index in lon_indices:
                    overlap_area = lat_old_cell_size * lon_old_cell_size  # Approximate overlap area
                    regridded_data_values[i, j] += burned_area_data[lat_index, lon_index].item() * overlap_area

    return regridded_data
# Earth surface area matrix for 0.25x0.25 degree cells (in m^2)
#earth_surface_area = np.full(old_shape, 510.1e12 / (old_nlon * old_nlat))  # Earth's surface area is ~510.1 million km^2

# Process each file
for filepath in nc_files_list:
    if gfed == '5':
        gfed_data = xr.open_dataset(filepath)
        # Regrid specified variables and calculate 'Nat'
        total_regridded = regrid_burned_area(gfed_data['Total'], target_lat, target_lon)
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
        #small_fire_area_regridded = regrid_data_xarray(small_fire_area, lat1, lon1)
        
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
