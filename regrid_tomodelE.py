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


# Range of years to process
years = range(2001, 2002)

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
    Perform conservative regridding from a high-resolution grid to a lower-resolution grid using vectorized operations.

    The method calculates overlaps between the source grid cells (high-resolution) and the target grid cells
    (lower-resolution). It then sums the values from the source grid cells based on these overlaps and assigns
    the total to the corresponding target grid cell to ensure conservation of the total sum.

    Special handling is applied for the polar regions, where the latitude ranges from -90 to -87 and from 87 to 90
    are treated as single cells with a size of 3 degrees.
    
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
    source_values = burned_area_data.values
    
    # Calculate the grid cell size for both grids
    lat_old_cell_size = 180./source_lat.size
    lon_old_cell_size = 360./source_lon.size
    lat_new_cell_size = 180./target_lat.size
    lon_new_cell_size = 360./target_lon.size 

    # Initialize the new regridded array
    regridded_data_values = np.zeros((target_lat.size, target_lon.size))

    # Convert target_lat and target_lon to NumPy arrays if they aren't already
    target_lat = np.asarray(target_lat)
    target_lon = np.asarray(target_lon)

    # Vectorized calculation of the min and max bounds for the target grid cells
    lat_min = target_lat - lat_new_cell_size / 2
    lat_max = target_lat + lat_new_cell_size / 2

    # Handle special cases for polar latitudes
    lat_min[target_lat <= -87] = -90.
    lat_max[target_lat >= 87] = 90.

    lon_min = target_lon - lon_new_cell_size / 2
    lon_max = target_lon + lon_new_cell_size / 2

    # Vectorized calculation of the min and max bounds for the source grid cells
    src_lat_min = source_lat - lat_old_cell_size / 2
    src_lat_max = source_lat + lat_old_cell_size / 2
    src_lon_min = source_lon - lon_old_cell_size / 2
    src_lon_max = source_lon + lon_old_cell_size / 2

    # Pre-calculate latitude overlaps once per target latitude
    lat_overlap = np.maximum(0, np.minimum(lat_max[:, np.newaxis], src_lat_max) - np.maximum(lat_min[:, np.newaxis], src_lat_min))

    # Calculate the overlaps and sum values for each target cell
    for j in range(target_lon.size):
        # Calculate the overlap in longitude
        lon_overlap = np.maximum(0, np.minimum(lon_max[j], src_lon_max) - np.maximum(lon_min[j], src_lon_min))

        # Compute the total overlap area for each target cell
        overlap_area = lat_overlap * lon_overlap

        # Sum the contributions from the source grid cells to the target grid cell
        regridded_data_values[:, j] = np.sum(source_values * overlap_area, axis=(1, 2))

    # Convert the regridded array back to an xarray DataArray with the correct coordinates
    regridded_data = xr.DataArray(regridded_data_values, coords=[target_lat, target_lon], dims=['lat', 'lon'])
    # Scale the regridded data to conserve total sum
    total_source_sum = np.sum(source_values)
    total_target_sum = np.sum(regridded_data_values)

    if total_target_sum != 0:
        regridded_data *= total_source_sum / total_target_sum

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
            peat_regridded = regrid_burned_area(gfed_data['Peat'], target_lat, target_lon)
            crop_regridded = regrid_burned_area(gfed_data['Crop'], target_lat, target_lon)
            defo_regridded = regrid_burned_area(gfed_data['Defo'], target_lat, target_lon)
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
                'Total': total_regridded
            })
        
    # Generate unique output file name
    output_filename = os.path.join(outputpath, f"regridded_{os.path.basename(filepath)}")
    # Save the regridded dataset to a new NetCDF file
    regridded_dataset.to_netcdf(output_filename)

    print(f"Regridded data saved to {output_filename}")
