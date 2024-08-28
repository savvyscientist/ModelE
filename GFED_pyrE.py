import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import subprocess
import numpy as np
import xarray as xr
import sys
import netCDF4 as nc
from netCDF4 import Dataset

################################################################

def load_config():
    """Load configuration settings for file paths, year range, and other parameters."""
    config = {
        'd': '/discover/nobackup/bkanzer/',
        'd2': '/discover/nobackup/bkanzer/nk_CCycle_E6obioF40',
        'd3': '/discover/nobackup/projects/giss/prod_input_files/emis/BBURN_ALT/20240517/BBURN_GFED_4s/monthly',
        'd31': '/discover/nobackup/projects/giss/prod_input_files/emis/BBURN_ALT/20240517/BBURN_GFED_4s/monthly/NAT',
        'd4': '/discover/nobackup/projects/giss_ana/users/nkiang/aDATA/Biomass/Spawn2020/Spawn_Ent',
        'd5': '/discover/nobackup/nkiang/DATA/Vegcov/V2HX2_EntGVSD_v1.1.2_16G_Spawn2020Sa_biomass_agb_2010_ann_pure.nc',
        'cmap': 'jet',
        'iyear': 2010,
        'fyear': 2010,
        'func_group': ["grass", "shrub", "tree"],
        'srang': ["009", "010"],
        'trang': ["002", "004", "006", "007", "008", "009", "010"],
        'grang': ["011", "012", "013", "014"]
    }
    return config

################################################################

def calculate_emis(vtype, BA, f, missing_val, nan_mat, srang, trang, grang, lats, lons, zero_mat, file_path):
    """Calculate emissions based on vegetation type and other factors."""
    d = Dataset(file_path, 'r') 
    
    looptype_mapping = {'grass': grang, 'shrub': srang, 'tree': trang}
    
    var = np.zeros((len(lats), len(lons)), dtype=float)
    vf_arr = np.zeros((len(lats), len(lons)), dtype=float)
    
    looptype = looptype_mapping[vtype]
    
    for t in looptype:
        vf = f.variables[f"ra001{t}"][:]
        vf = np.where(vf < 0., zero_mat, vf)
        vf_arr += vf

        totCM = f.variables[f"ra017{t}"][:]
        totCM = np.where(totCM < 0., zero_mat, totCM)

        soilC = f.variables[f"ra024{t}"][:]
        soilC = np.where(soilC < 0., zero_mat, soilC)

        diff = (totCM - soilC)
        
        var += (vf * (diff))
    
    emis_type = np.divide(BA * var, vf_arr, where=vf_arr != 0, out=np.full_like(vf_arr, 0.))
    
    return emis_type

################################################################

def calculate_spawnemis(vtype, BA, zero_mat, d, f, srang, trang, grang, lats, lons, file_path):
    """Calculate spawn emissions based on vegetation type and other factors."""
    d = Dataset(file_path, 'r')
    
    looptype_mapping = {'grass': grang, 'shrub': srang, 'tree': trang}
    looptype = looptype_mapping[vtype]
    
    vf_arr = np.zeros((len(lats), len(lons)), dtype=float)
    mult_arr = np.zeros((len(lats), len(lons)), dtype=float)

    for t in looptype:
        vf = f.variables[f"ra001{t}"][:]
        vf = np.where(vf < 0., zero_mat, vf)
        vf_arr += vf
        
        SpawnCM_ABV = d.variables['biomass'][t,:,:]
        mult = vf * SpawnCM_ABV
        mult_arr += mult

    emis_type = np.divide((BA * mult_arr), vf_arr, where=vf_arr != 0, out=np.full_like(vf_arr, 0.))

    return emis_type

################################################################

def define_subplot(ax, data, lons, lats, cmap, cborientation, fraction, pad, labelpad, fontsize, title, glob, clabel, masx, is_diff=False):
    """Define the properties of a subplot with optional difference normalization."""
    ax.coastlines(color='black')
    ax.add_feature(cfeature.LAND, edgecolor='gray')
    ax.add_feature(cfeature.OCEAN, facecolor='white', edgecolor='none', zorder=1)

    ax.set_title(title, fontsize=10, pad=25)
    props = dict(boxstyle="round", facecolor='lightgray', alpha=0.5)
    ax.text(0.5, 1.07, f"Global Total: {glob}", ha="center", va="center", transform=ax.transAxes, bbox=props, fontsize=10)

    if is_diff:
        data_min, data_max = data.min(), data.max()
        if data_min == data_max:
            norm = mcolors.Normalize(vmin=data_min - 1, vmax=data_max + 1)
        else:
            abs_max = max(abs(0.25 * data_min), abs(0.25 * data_max))
            norm = mcolors.Normalize(vmin=-abs_max, vmax=abs_max)
    else:
        norm = None

    p = ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm, vmin=0 if not is_diff else None, vmax=masx if not is_diff else None)
    cbar = fig.colorbar(p, ax=ax, orientation=cborientation, fraction=fraction, pad=pad)
    cbar.set_label(clabel, labelpad=labelpad, fontsize=fontsize)

    return ax

################################################################

def process_data():
import matplotlib.pyplot as plt

def process_data():
    """Main data processing function that integrates various components."""
    config = load_config()
    
    emis_G = np.zeros((90, 144), dtype=float)
    emis_S = np.zeros((90, 144), dtype=float)
    emis_T = np.zeros((90, 144), dtype=float)

    for year in range(config['iyear'], config['fyear'] + 1):
        filename = f"ANN{year}.aijnk_CCycle_E6obioF40.nc"
        filepath = os.path.join(config['d'], filename)
    
        if os.path.exists(filepath):
            with nc.Dataset(filepath) as f:
                lats = f.variables['lat'][:]
                lons = f.variables['lon'][:]
                
                # For each functional group, calculate the emissions using the main dataset
                for vtype in config['func_group']:
                    if vtype == 'grass':
                        emis_G += calculate_emis(vtype, BA=None, f=f, missing_val=None, nan_mat=None, 
                                                 srang=config['srang'], trang=config['trang'], 
                                                 grang=config['grang'], lats=lats, lons=lons, 
                                                 zero_mat=np.zeros_like(lats), file_path=filepath)
                    elif vtype == 'shrub':
                        emis_S += calculate_emis(vtype, BA=None, f=f, missing_val=None, nan_mat=None, 
                                                 srang=config['srang'], trang=config['trang'], 
                                                 grang=config['grang'], lats=lats, lons=lons, 
                                                 zero_mat=np.zeros_like(lats), file_path=filepath)
                    elif vtype == 'tree':
                        emis_T += calculate_emis(vtype, BA=None, f=f, missing_val=None, nan_mat=None, 
                                                 srang=config['srang'], trang=config['trang'], 
                                                 grang=config['grang'], lats=lats, lons=lons, 
                                                 zero_mat=np.zeros_like(lats), file_path=filepath)
                
                # Calculate emissions using the Spawn dataset and combine them with the main dataset emissions
                spawn_filepath = config['d5']  # Path to the Spawn dataset
                emis_G += calculate_spawnemis(vtype='grass', BA=None, zero_mat=np.zeros_like(lats), 
                                              d=f, f=f, srang=config['srang'], 
                                              trang=config['trang'], grang=config['grang'], 
                                              lats=lats, lons=lons, file_path=spawn_filepath)
                
                emis_S += calculate_spawnemis(vtype='shrub', BA=None, zero_mat=np.zeros_like(lats), 
                                              d=f, f=f, srang=config['srang'], 
                                              trang=config['trang'], grang=config['grang'], 
                                              lats=lats, lons=lons, file_path=spawn_filepath)
                
                emis_T += calculate_spawnemis(vtype='tree', BA=None, zero_mat=np.zeros_like(lats), 
                                              d=f, f=f, srang=config['srang'], 
                                              trang=config['trang'], grang=config['grang'], 
                                              lats=lats, lons=lons, file_path=spawn_filepath)

        else:
            print(f"File {filepath} not found. Skipping.")
    
    # Step 3: Visualize the Results
    fig, ax = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    
    define_subplot(ax[0], emis_G, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Grass Emissions', 
                   glob=None, clabel='Emissions', masx=emis_G.max())
    
    define_subplot(ax[1], emis_S, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Shrub Emissions', 
                   glob=None, clabel='Emissions', masx=emis_S.max())
    
    define_subplot(ax[2], emis_T, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Tree Emissions', 
                   glob=None, clabel='Emissions', masx=emis_T.max())

    plt.show()



################################################################

# Run the main process
if __name__ == "__main__":
    process_data()

