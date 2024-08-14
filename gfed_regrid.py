import h5py
import numpy as np
from scipy import interpolate
import xarray as xr
from scipy.ndimage import zoom
from netCDF4 import Dataset
import sys
import os
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

gfed='4s'
if gfed == '4s':
    filepath = '/discover/nobackup/projects/giss_ana/users/kmezuman/gfed4s/updated'
elif gfed == '5':
    filepath = '/discover/nobackup/projects/giss_ana/users/kmezuman/GFED5/BA'
else:
    stop,'code not compatible with GFED f"{gfed}"'

latpath = '/discover/nobackup/kemzuman/ANN2005.aijnk_CCycle_E6obioF40.nc'

new_shape = (90, 144)
old_shape = (180, 360)
output_base = filepath+'/regrid'

filename_group = []
filepath_group = []

dataset = 'burned_area'

years = range(1997,2020)

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
        stop,'code not compatible with GFED f"{gfed}"'
    filename_group.append(filename)


for m in filename_group:
    filepath_place = os.path.join(filepath, m)
    filepath_group.append(filepath_place)

with Dataset(latpath, 'r') as d:
    if 'lat' in d.variables and 'lon' in d.variables:
        lat1 = d.variables['lat'][:]
        lon1 = d.variables['lon'][:]
    else:
        raise KeyError("Variables 'lat' and 'lon' not found in the NetCDF file.")

for idx in filepath_group:
    try:
        with nc.Dataset(idx, 'r') as file:
            peat = crop = defo = total = None

            if 'Peat' in file.variables:
                peat = file['Peat'][:]
            if 'Crop' in file.variables:
                crop = file['Crop'][:]
            if 'Defo' in file.variables:
                defo = file['Defo'][:]
            if 'Total' in file.variables:
                total = file['Total'][:]

            peat = np.ma.filled(peat, 0)
            crop = np.ma.filled(crop, 0)
            defo = np.ma.filled(defo, 0)
            total = np.ma.filled(total, 0)

            peat = peat[0, :, :]
            crop = crop[0, :, :]
            defo = defo[0, :, :]
            total = total[0, :, :]

            burnedarea = total - (peat + defo + crop)

            ds = xr.Dataset(
                data_vars={
                    'burnedarea': (['lat', 'lon'], burnedarea)
                },
                coords={
                    'lat': (['lat'], np.linspace(-90, 90, burnedarea.shape[0])),
                    'lon': (['lon'], np.linspace(-180, 180, burnedarea.shape[1]))
                }
            )
            
            target_lat = np.linspace(-90, 90, 90)
            target_lon = np.linspace(-180, 180, 144)

            region_data_regridded = ds['burnedarea'].interp(lat=target_lat, lon=target_lon, method='nearest')
            resized_data = region_data_regridded.values

            #output_file = f"{output_base}{idx + 1}.nc"
            output_file = f"{output_base}{os.path.basename(idx).split('.')[0]}.nc"
            if os.path.exists(output_file):
                with Dataset(output_file, 'w', format='NETCDF4') as outfile:
                    outfile.createDimension('lat', len(target_lat))
                    outfile.createDimension('lon', len(target_lon))
                    resized_BA = outfile.createVariable('burnedarea', 'f4', ('lat', 'lon'))
                    resized_BA[:] = resized_data

                    lat_var = outfile.createVariable('lat', 'f4', ('lat',))
                    lon_var = outfile.createVariable('lon', 'f4', ('lon',))
                    lat_var[:] = target_lat
                    lon_var[:] = target_lon

                    resized_BA.units = 'm$^2$'
                    resized_BA.long_name = 'Burned Area'
                    lat_var.units = 'degrees_north'
                    lon_var.units = 'degrees_east'

                    from pathlib import Path
                    file_name = Path(output_file).name

                    print(f"Data successfully written to NetCDF file titled {file_name}.")


    except KeyError as e:
        print(f"KeyError: {e} not found in file")
    except ValueError as e:
        print(f"ValueError: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")    
