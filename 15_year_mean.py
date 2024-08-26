import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import subprocess
import numpy as np
import xarray as xr
import sys
import netCDF4 as nc
import matplotlib.colors as mcolors
################################################################

def calculate_emis(vtype, BA, f, missing_val, nan_mat):
    srang = ["009", "010"]
    trang = ["001", "002", "003", "004", "005", "006", "007", "008", "009", "010"]
    grang = ["011", "012", "013", "014"]
    looptype_mapping = {'grass': grang, 'shrub': srang, 'tree': trang}
    var = np.zeros((len(lats), len(lons)),dtype=float)
    vf_arr = np.zeros((len(lats), len(lons)),dtype=float)
    looptype = looptype_mapping[vtype]

    for t in looptype:
        vf = f.variables[f"ra001{t}"][:]
        #vf = np.where(vf <= missing_val, nan_mat, vf)
        vf_arr += vf
        ###print(len(vf))
        #sys.exit()
        totCM = f.variables[f"ra017{t}"][:]
        #totCM = np.where(totCM <= missing_val, nan_mat, totCM)

        soilC = f.variables[f"ra024{t}"][:]
        #soilC = np.where(soilC <= missing_val, nan_mat, soilC)
        var += (vf * (totCM - soilC))
    if np.all(vf_arr==0):
        pass
    else:
        emis_type = BA * var / vf_arr
    #print(emis_type)
    return emis_type

################################################################
#fraction=0.064, pad=0.15, labelpad=-45, fontsize=10
def define_subplot(ax, data, lons, lats, cmap, cborientation, fraction, pad, labelpad, fontsize, title, clabel):
    ax.coastlines(color='black')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.LAND, edgecolor='gray')

    ax.set_title(title)

    p = ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(), cmap=cmap)
    cbar = fig.colorbar(p, ax=ax, orientation=cborientation, fraction=fraction, pad=pad)
    cbar.set_label(clabel, labelpad=labelpad, fontsize=fontsize)
    return ax
#Offline kgCyr^{-1}
################################################################

def define_diff_subplot(ax, data, lons, lats, cmap, cborientation, fraction, pad, labelpad, fontsize, title, clabel):
    ax.coastlines(color='black')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.LAND, edgecolor='gray')

    ax.set_title(title)

    data_min, data_max = data.min(), data.max()
    norm = mcolors.Normalize(vmin=-abs(data_max), vmax=abs(data_max))

    p = ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(p, ax=ax, orientation=cborientation, fraction=fraction, pad=pad)
    cbar.set_label(clabel, labelpad=labelpad, fontsize=fontsize)
    return ax

################################################################ 

d = '/discover/nobackup/bkanzer/'
d2 = '/discover/nobackup/bkanzer/nk_CCycle_E6obioF40'
d3 = '/discover/nobackup/projects/giss/prod_input_files/emis/BBURN_ALT/20240517/BBURN_GFED_4s/monthly'

cmap = 'jet' #cividis, Spectral, jet
iyear = 2000
fyear = 2015
func_group = ['grass','shrub','tree']

################################################################

emis_G = 0
emis_S = 0
emis_T = 0

emis15_sum = np.zeros((90, 144), dtype=float)
cO2v15_sum = np.zeros((90, 144), dtype=float)
GFED15_sum = np.zeros((90, 144), dtype=float)
mean15_emis = np.zeros((90, 144), dtype=float)

for object in range(1):
    for year in range(iyear, fyear + 1):
        filename = f"ANN{year}.aijnk_CCycle_E6obioF40.nc"
        filepath = os.path.join(d, filename)
        if os.path.exists(filepath):
            with nc.Dataset(filepath) as f:
                emis15_sum = 0
                lats = f.variables['lat'][:]
                lons = f.variables['lon'][:]
                missing_val = -1E-30
                nan_mat = np.full((len(lats), len(lons)), np.nan)
                zero_mat = np.zeros((len(lats), len(lons)), dtype=float)
                for vtype in func_group:
                    if vtype == "grass":
                        BA_G = f.variables['BA_grass'][:]
                        #BA_G = np.where(BA_G <= missing_val, nan_mat, BA_G)
                        emis_G = calculate_emis(vtype, BA_G, f, missing_val, nan_mat)
                    elif vtype == "shrub":
                        BA_S = f.variables['BA_shrub'][:]
                        #BA_S = np.where(BA_S <= missing_val, nan_mat, BA_S)
                        emis_S = calculate_emis(vtype, BA_S, f, missing_val, nan_mat)
                    if vtype == "tree":
                        BA_T = f.variables['BA_tree'][:]
                        #BA_T = np.where(BA_T <= missing_val, nan_mat, BA_T)
                        emis_T = calculate_emis(vtype, BA_T, f, missing_val, nan_mat)
                    else:
                        break
                emis_tot = emis_G + emis_S + emis_T
                emis_min = emis_tot[:].min()
                emis_max = emis_tot[:].max()
            emis15_sum += emis_tot
        mean15_emis += (emis15_sum / 15.)

    filename2 = f"ANN{year}.taijnk_CCycle_E6obioF40.nc"
    filepath2 = os.path.join(d2, filename2)

    if os.path.exists(filepath2):
        cO2v_data = []
        with nc.Dataset(filepath2) as f2:
            lats2 = f2.variables['lat'][:]
            lons2 = f2.variables['lon'][:]
            cO2v = f2.variables['CO2n_pyrE_src'][:]
            #cO2v = np.where(cO2v <= 0., nan_mat, cO2v)

            vmin2 = cO2v[:].min()
            vmax2 = cO2v[:].max()
            cO2v *= 1.0E-12
            axyp = f2.variables['axyp'][:]

            cO2v *= axyp
            cO2v *= (60.*60.*24.*365.25)

            cO2v15_sum += cO2v

        cO2v15_mean = cO2v15_sum / 12.
        
        ds15_mean = mean15_emis - cO2v15_mean
        vmin3 = ds15_mean[:].min()
        vmax3 = ds15_mean[:].max()
    
    filename3 = f"{year}.nc"
    filepath3 = os.path.join(d3, filename3)
    if os.path.exists(filepath3):
        missing_val2 = 0.
        ann_sum = np.zeros((len(lats2), len(lons2)), dtype=float)
        with nc.Dataset(filepath3) as f3:
            for year in range(1): #previously did range(15):
                for k in range(12):

                    GFED_data = f3.variables['CO2n'][k, :, :]
                    GFED_CO2 = GFED_data.reshape(90, 144)

                    GFED_CO2 = np.where(GFED_CO2 <= 0., zero_mat, GFED_CO2)
                    GFED_CO2  *= axyp

                    GFED_CO2  *= (60.*60.*24.*365.25) #If I do the sum I don't need this conversion... ?
                ann_sum += GFED_CO2  
            GFann_mean = ann_sum / 12.

            GFED15_sum += GFann_mean
            GFED15_mean = GFED15_sum / 15.
            diff_srdata = GFED15_mean - cO2v15_mean
            diff_rsdata = GFED15_mean - mean15_emis

################################################################

tmin = min(0, emis_min, vmin2, vmin3)
tmax = max(500, emis_max, vmax2, vmax3)

################################################################

frac = 0.05
pad = 0.025
lp = 5
fs = 10
fig, axs = plt.subplots(3, 2, figsize=(50, 20), subplot_kw={'projection': ccrs.PlateCarree()})  # Use GeoAxes for all
#define_subplot(axs[0, 0], emis_tot, lons, lats, cmap, "horizontal", 0.064, 0.15, -45, 10, f"Offline Emissions {year}", 'Emissions [$kgCyr^{-1}$]')
define_subplot(axs[0, 0], mean15_emis, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Mean Offline Emissions {iyear} - {fyear}", 'Emissions [$kgCyr^{-1}$]')
define_subplot(axs[1, 0], cO2v15_mean, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Mean Online Emissions {iyear} - {fyear}", 'Emissions [$kgyr^{-1}$]')
define_subplot(axs[0, 1], GFED15_mean, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"GFED Observed Mean Emissions {iyear} - {fyear}", 'Emissions [$kgyr^{-1}$]')
define_diff_subplot(axs[2, 0], ds15_mean, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Offline and Online", 'Emissions [$kgCyr^{-1}-kgyr^{-1}$]')
define_diff_subplot(axs[1, 1], diff_srdata, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of GFED Observed Emissions and ModelE Online", 'Emissions [$kgCyr^{-1}$]')
define_diff_subplot(axs[2, 1], diff_rsdata, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of GFED Observed Emissions and Offline Calculation", 'Emissions [$kgCyr^{-1}$]')

#fig.suptitle(' 15 Year Mean', fontsize=30)
plt.tight_layout(pad=3.0, w_pad=2.0, h_pad=7.0)
plt.show()
