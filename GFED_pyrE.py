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
import re

################################################################

def load_config():
    """Load configuration settings for file paths, year range, and other parameters."""
    config = {
        'dir_sim': '/discover/nobackup/kmezuman/nk_CCycle_E6obioF40',
        'dir_obs_emis': '/discover/nobackup/projects/giss/prod_input_files/emis/BBURN_ALT/20240517/BBURN_GFED_4s/monthly/NAT',
        'dir_obs_bio': '/discover/nobackup/nkiang/DATA/Vegcov/V2HX2_EntGVSD_v1.1.2_16G_Spawn2020Sa_biomass_agb_2010_ann_pure.nc',
        'nlat': 90,
        'nlon': 144,
        'cmap': 'jet',
        'iyear': 2010,
        'fyear': 2010,
        'srang': ["009", "010"],
        'trang': ["002", "004", "006", "007", "008", "009", "010"],
        'grang': ["011", "012", "013", "014"]
    }
    return config

################################################################

def extract_scaling_factor(units):
    """
    Extracts the scaling factor from the units string and converts it to a number.
    For example, '10^-3 kg/m^2/s' will return 0.001 and 'kg/m^2/s'.
    """
    match = re.match(r"^(10\^(-?\d+)|[-+]?\d*\.?\d+([eE][-+]?\d+)?)\s*(.*)$", units)
    if match:
        if match.group(1).startswith("10^"):
            scaling_factor = float(10) ** float(match.group(2))
        else:
            scaling_factor = float(match.group(1))
        unit = match.group(4)
        return scaling_factor, unit
    return 1.0, units  # Default scaling factor is 1 if not specified

################################################################

def calculate_emis(vrang, BA, f, missing_val, nan_mat, lats, lons, zero_mat):
    """Calculate emissions based on vegetation type and other factors."""
    
    var = np.zeros((len(lats), len(lons)), dtype=float)
    vf_arr = np.zeros((len(lats), len(lons)), dtype=float)
    
    for t in vrang:
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

def calculate_spawnemis(vrang, BA, file_path, f, lats, lons, zero_mat):
    """Calculate spawn emissions based on vegetation type and other factors."""
    d = Dataset(file_path, 'r') 
    vf_arr = np.zeros((len(lats), len(lons)), dtype=float)
    mult_arr = np.zeros((len(lats), len(lons)), dtype=float)

    for t in vrang:
        vf = f.variables[f"ra001{t}"][:]
        vf = np.where(vf < 0., zero_mat, vf)
        vf_arr += vf
        
        SpawnCM_ABV = d.variables['biomass_agb'][t,:,:]
        mult = vf * SpawnCM_ABV
        mult_arr += mult

    emis_type = np.divide((BA * mult_arr), vf_arr, where=vf_arr != 0, out=np.full_like(vf_arr, 0.))

    return emis_type

################################################################

def define_subplot(fig, ax, data, lons, lats, cmap, cborientation, fraction, pad, labelpad, fontsize, title, glob, clabel, masx, is_diff=False):
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
    zero_mat = np.zeros((config['nlat'], config['nlon']), dtype=float)
    spawn_filepath = config['dir_obs_bio']
    

    # Calculate emissions from 1. Ent Biomass and  2. Spawn biomass
    # using pyrE BA and assumption of 100% combustion completeness
    for year in range(config['iyear'], config['fyear'] + 1):
        veg_filename = f"ANN{year}.aijnk_CCycle_E6obioF40.nc"
        filepath = os.path.join(config['dir_sim'], veg_filename)
        zero_mat = np.zeros((config['nlat'], config['nlon']), dtype=float)
    
        if os.path.exists(filepath):
            # Initialize arrays for emissions from different datasets
            emis_G = np.zeros((config['nlat'], config['nlon']), dtype=float)
            emis_S = np.zeros((config['nlat'], config['nlon']), dtype=float)
            emis_T = np.zeros((config['nlat'], config['nlon']), dtype=float)
            emis_sp_G = np.zeros((config['nlat'], config['nlon']), dtype=float)
            emis_sp_S = np.zeros((config['nlat'], config['nlon']), dtype=float)
            emis_sp_T = np.zeros((config['nlat'], config['nlon']), dtype=float)
            with nc.Dataset(filepath) as f:
                lats = f.variables['lat'][:]
                lons = f.variables['lon'][:]
                axyp = f.variables['axyp'][:]
                
                # Calculate emissions using Ent biomass and pyrE BA
                # and calculate emissions using Spawn biomass and pyrE BA
                BA_grass = f.variables['BA_grass']

                emis_G = calculate_emis(vrang=config['grang'], BA=BA_grass, f=f, missing_val=None, nan_mat=None, 
                                         lats=lats, lons=lons, zero_mat=zero_mat)

                emis_sp_G = calculate_spawnemis(vrang=config['grang'], BA=BA_grass, file_path=spawn_filepath,
                                                 f=f, lats=lats, lons=lons, zero_mat=zero_mat)

                BA_shrub = f.variables['BA_shrub']

                emis_S = calculate_emis(vrang=config['srang'], BA=BA_shrub, f=f, missing_val=None, nan_mat=None, 
                                         lats=lats, lons=lons, zero_mat=zero_mat)

                emis_sp_S = calculate_spawnemis(vrang=config['srang'], BA=BA_shrub, file_path=spawn_filepath,
                                                 f=f, lats=lats, lons=lons, zero_mat=zero_mat)

                BA_tree = f.variables['BA_tree']

                emis_T = calculate_emis(vrang=config['trang'], BA=BA_tree, f=f, missing_val=None, nan_mat=None, 
                                         lats=lats, lons=lons, zero_mat=zero_mat)

                emis_sp_S = calculate_spawnemis(vrang=config['trang'], BA=BA_tree, file_path=spawn_filepath,
                                                 f=f, lats=lats, lons=lons, zero_mat=zero_mat)
                
                emis_tot = emis_G + emis_S + emis_T   # Units kgCm^(-2)yr^(-1)
                tot_emis = (np.nansum(emis_tot))*60.*60.*24.*365.
                tot_emis = format(tot_emis, '.3e')

                emis_sp_tot = emis_sp_G + emis_sp_S + emis_sp_T   # Units kgCm^(-2)yr^(-1)
                tot_sp_emis = (np.nansum(emis_sp_tot))*60.*60.*24.*365.
                tot_sp_emis = format(tot_sp_emis, '.3e')
        else:
            print(f"File {filepath} not found. Skipping.")

        # Read pyrE CO2n emissions
        tracer_filename = f"ANN{year}.taijnk_CCycle_E6obioF40.nc"
        filepath = os.path.join(config['dir_sim'], tracer_filename)
        if os.path.exists(filepath):
            with nc.Dataset(filepath) as f:
                emis_pyrE = f.variables['CO2n_pyrE_src'][:]
                units = f.variables['CO2n_pyrE_src'].units
                scaling_factor, unit = extract_scaling_factor(units)
                emis_pyrE *= float(scaling_factor)
        else:
            print(f"File {filepath} not found. Skipping.")

        # Calculate GFED4s annual mean emissions for comparison
        obs_filepath = os.path.join(config['dir_obs_emis'], f"{year}.nc")  # Example: Use actual filenames here
    
        if os.path.exists(obs_filepath):
            ann_sum = np.zeros((config['nlat'], config['nlon']), dtype=float)
            with nc.Dataset(obs_filepath) as f_obs:
                for k in range(12):
                    GFED_data = f_obs.variables['CO2n'][k, :, :]  
                    GFED_CO2 = GFED_data.reshape(config['nlat'], config['nlon'])
                    GFED_CO2  *= axyp
                    GFED_CO2  *= (60.*60.*24.*365.) 
                    GFED_CO2 = np.where(GFED_CO2 <= 0., zero_mat, GFED_CO2)
                    ann_sum += GFED_CO2
                totGFED = np.nansum(ann_sum)
                totGFED = format(totGFED, '3e')
                ann_mean = ann_sum / 12.

        else:
            print(f"File {obs_filepath} not found. Skipping.")

    #Once the GFED4s and GFED5 BA is regridded emissions can be calculated with both the Ent and Spawn biomass
    
    # Step 3: Visualize the Results
    fig, ax = plt.subplots(2, 2, figsize=(18, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Plot emissions from pyrE BA and Ent Biomass 
    define_subplot(fig, ax[0, 0], emis_tot, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Offline Emissions (pyrE BA and Ent Biomass)', 
                   glob=None, clabel='Emissions', masx=0.7*emis_tot.max())
    
    # Plot emissions from the main dataset and Spawn dataset
    define_subplot(fig, ax[0, 1], emis_sp_tot, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='Offline Emissions (pyrE BA and Spawn Biomass)', 
                   glob=None, clabel='Emissions', masx=0.7*emis_sp_tot.max())
    
    # Plot emissions from pyrE 
    define_subplot(fig, ax[1, 0], emis_pyrE, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='pyrE old emissions', 
                   glob=None, clabel='Emissions', masx=0.7*emis_pyrE.max())

    # Plot emissions from GFED4s
    define_subplot(fig, ax[1, 1], ann_mean, lons, lats, cmap=config['cmap'], 
                   cborientation='horizontal', fraction=0.05, pad=0.05, 
                   labelpad=5, fontsize=10, title='GFED4s emissions', 
                   glob=None, clabel='Emissions', masx=0.7*ann_mean.max())
    
    plt.tight_layout()
    plt.show()

################################################################

# Run the main process
if __name__ == "__main__":
    process_data()

