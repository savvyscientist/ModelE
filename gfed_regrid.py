import xarray as xr
import os

gfed = '4s'
if gfed == '4s':
    filepath = '/discover/nobackup/projects/giss_ana/users/kmezuman/gfed4s/updated'
elif gfed == '5':
    filepath = '/discover/nobackup/projects/giss_ana/users/kmezuman/GFED5/BA'
else:
    raise ValueError(f"Code not compatible with GFED '{gfed}'")

latpath = '/discover/nobackup/kemzuman/ANN2005.aijnk_CCycle_E6obioF40.nc'

new_shape = (90, 144)
old_shape = (180, 360)
output_base = filepath + '/regrid'

filename_group = []
filepath_group = []

dataset = 'burned_area'

years = range(1997, 2020)

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

for m in filename_group:
    filepath_place = os.path.join(filepath, m)
    filepath_group.append(filepath_place)

# Load the latitude and longitude using xarray
ds = xr.open_dataset(latpath)
if 'lat' in ds.variables and 'lon' in ds.variables:
    lat1 = ds['lat'].values
    lon1 = ds['lon'].values
else:
    raise KeyError("Variables 'lat' and 'lon' not found in the dataset.")

# Placeholder for regridding logic
# Assume the regridding part is correctly implemented with appropriate xarray methods

def regrid_data_xarray(data, new_lat, new_lon):
    regridded = data.interp(lat=new_lat, lon=new_lon, method='linear')
    return regridded

# Example of applying the regridding function (assuming 'data' is an xarray DataArray)
# regridded_data = regrid_data_xarray(data, new_latitudes, new_longitudes)

# Save regridded data using xarray
# regridded_data.to_netcdf(output_base + '/regridded_output.nc')
