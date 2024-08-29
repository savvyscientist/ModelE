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
        'd': '/discover/nobackup/kmezuman/',
        'd2': '/discover/nobackup/kmezuman/nk_CCycle_E6obioF40',
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
    """Main data processing function that integrates various components."""
    config = load_config()
    
    # Initialize arrays for emissions from different datasets
    emis_G = np.zeros((90, 144), dtype=float)
    emis_S = np.zeros((90, 144), dtype=float)
    emis_T = np.zeros((90, 144), dtype=float)

    # Process emissions from the main dataset and Spawn dataset
    for year in range(config['iyear'], config['fyear'] + 1):
        filename = f"ANN{year}.aijnk_CCycle_E6obioF40.nc"
        filepath = os.path.join(config['d'], filename)
    
        if os.path.exists(filepath):
            with nc.Dataset(filepath) as f:
                lats = f.variables['lat'][:]
                lons = f.variables['lon'][:]
                
                # Calculate emissions using the main dataset
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
                
                # Calculate emissions using the Spawn dataset (d5)
                spawn_filepath = config['d5']
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
    
    # Processing emissions from the d2 dataset (e.g., additional processing or transformation steps)
    d2_filepath = os.path.join(config['d2'], 'some_specific_file.nc')  # Example: Use actual filenames here
    
    if os.path.exists(d2_filepath):
        with nc.Dataset(d2_filepath) as f_d2:
            lats = f_d2.variables['lat'][:]
            lons = f_d2.variables['lon'][:]
            
            # Example of processing specific to d2 dataset (e.g., filtering or special handling)
            emis_G_d2 = np.zeros((len(lats), len(lons)), dtype=float)
            emis_S_d2 = np.zeros((len(lats), len(lons)), dtype=float)
            emis_T_d2 = np.zeros((len(lats), len(lons)), dtype=float)
            
            for vtype in config['func_group']:
                emis_G_d2 += calculate_emis(vtype, BA=None, f=f_d2, missing_val=None, nan_mat=None, 
                                            srang=config['srang'], trang=config['trang'], 
                                            grang=config['grang'], lats=lats, lons=lons, 
                                            zero_mat=np.zeros_like(lats), file_path=d2_filepath)

    else:
        print(f"File {d2_filepath} not found. Skipping.")
    
    # Additional processing for d3 dataset
    d3_filepath = os.path.join(config['d3'], 'some_specific_file.nc')  # Example: Use actual filenames here
    
    if os.path.exists(d3_filepath):
        with nc.Dataset(d3_filepath) as f_d3:
            lats = f_d3.variables['lat'][:]
            lons = f_d3.variables['lon'][:]
            
            emis_G_d3 = np.zeros((len(lats), len(lons)), dtype=float)
            emis_S_d3 = np.zeros((len(lats), len(lons)), dtype=float)
            emis_T_d3 = np.zeros((len(lats), len(lons)), dtype=float)
            
            for vtype in config['func_group']:
                emis_G_d3 += calculate_emis(vtype, BA=None, f=f_d3, missing_val=None, nan_mat=None, 
                                            srang=config['srang'], trang=config['trang'], 
                                            grang=config['grang'], lats=lats, lons=lons, 
                                            zero_mat=np.zeros_like(lats), file_path=d3_filepath)

    else:
        print(f"File {d3_filepath} not found. Skipping.")
    
    # Processing emissions using d4 dataset (additional biomass or vegetation data)
    d4_filepath = os.path.join(config['d4'], 'some_other_specific_file.nc')  # Adjust as needed
    
    if os.path.exists(d4_filepath):
        with nc.Dataset(d4_filepath) as f_d4:
            lats = f_d4.variables['lat'][:]
            lons = f_d4.variables['lon'][:]
            
            emis_G_d4 = np.zeros((len(lats), len(lons)), dtype=float)
            emis_S_d4 = np.zeros((len(lats), len(lons)), dtype=float)
            emis_T_d4 = np.zeros((len(lats), len(lons)), dtype=float)
            
            for vtype in config['func_group']:
                emis_G_d4 += calculate_emis(vtype, BA=None, f=f_d4, missing_val=None, nan_mat=None, 
                                            srang=config['srang'], trang=config['trang'], 
                                            grang=config['grang'], lats=lats, lons=lons, 
                                            zero_mat=np.zeros_like(lats), file_path=d4_filepath)

    else:
        print(f"File {d4_filepath} not found. Skipping.")

    # Step 3: Visualize the Results
    fig, ax = plt.subplots(4, 3, figsize=(18, 16), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Plot emissions from the main dataset and Spawn dataset
    define_subplot(ax[0, 0], emis_G, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Grass Emissions (Main + Spawn)', 
                   glob=None, clabel='Emissions', masx=emis_G.max())
    
    define_subplot(ax[0, 1], emis_S, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Shrub Emissions (Main + Spawn)', 
                   glob=None, clabel='Emissions', masx=emis_S.max())
    
    define_subplot(ax[0, 2], emis_T, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Tree Emissions (Main + Spawn)', 
                   glob=None, clabel='Emissions', masx=emis_T.max())

    # Plot emissions from the d2 dataset
    define_subplot(ax[1, 0], emis_G_d2, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Grass Emissions (d2 Dataset)', 
                   glob=None, clabel='Emissions', masx=emis_G_d2.max())
    
    define_subplot(ax[1, 1], emis_S_d2, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Shrub Emissions (d2 Dataset)', 
                   glob=None, clabel='Emissions', masx=emis_S_d2.max())
    
    define_subplot(ax[1, 2], emis_T_d2, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Tree Emissions (d2 Dataset)', 
                   glob=None, clabel='Emissions', masx=emis_T_d2.max())

    # Plot emissions from the d3 dataset (BBURN_ALT dataset)
    define_subplot(ax[2, 0], emis_G_d3, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Grass Emissions (d3 Dataset)', 
                   glob=None, clabel='Emissions', masx=emis_G_d3.max())
    
    define_subplot(ax[2, 1], emis_S_d3, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Shrub Emissions (d3 Dataset)', 
                   glob=None, clabel='Emissions', masx=emis_S_d3.max())
    
    define_subplot(ax[2, 2], emis_T_d3, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Tree Emissions (d3 Dataset)', 
                   glob=None, clabel='Emissions', masx=emis_T_d3.max())
    
    # Plot emissions from the d4 dataset
    define_subplot(ax[3, 0], emis_G_d4, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Grass Emissions (d4 Dataset)', 
                   glob=None, clabel='Emissions', masx=emis_G_d4.max())
    
    define_subplot(ax[3, 1], emis_S_d4, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Shrub Emissions (d4 Dataset)', 
                   glob=None, clabel='Emissions', masx=emis_S_d4.max())
    
    define_subplot(ax[3, 2], emis_T_d4, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Tree Emissions (d4 Dataset)', 
                   glob=None, clabel='Emissions', masx=emis_T_d4.max())

    plt.tight_layout()
    plt.show()

################################################################

# Run the main process
if __name__ == "__main__":
    process_data()

