##At some point I stopped trying to make the code efficient. I apologize for that. It is a trend throughout a lot of my code.
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
#np.set_printoptions(threshold=sys.maxsize)  
################################################################

def calculate_emis(vtype, BA, f, missing_val, nan_mat, srang, trang, grang):
    d = Dataset(r'/discover/nobackup/bkanzer/ANN2003.aijnk_CCycle_E6obioF40.nc', 'r') 
    looptype_mapping = {'grass': grang, 'shrub': srang, 'tree': trang}
    var = np.zeros((len(lats), len(lons)),dtype=float)
    vf_arr = np.zeros((len(lats), len(lons)),dtype=float)
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

def calculate_spawnemis(vtype, BA, zero_mat, d, f, srang, trang, grang):
    d = Dataset(r'/discover/nobackup/projects/giss_ana/users/nkiang/aDATA/Biomass/Spawn2020/Spawn_Ent/V2HX2_EntGVSD_v1.1.2_16G_ann_pure_Spawn_AGB_2010.nc', 'r')
    looptype_mapping = {'grass': grang, 'shrub': srang, 'tree': trang}                                               
    looptype = looptype_mapping[vtype]
    vf_arr = np.zeros((len(lats), len(lons)),dtype=float)
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

def define_subplot(ax, data, lons, lats, cmap, cborientation, fraction, pad, labelpad, fontsize, title, glob, clabel, masx):
    ax.coastlines(color='black')
    ax.add_feature(cfeature.LAND, edgecolor='gray')

    ocean_feature = cfeature.OCEAN
    ax.add_feature(ocean_feature, facecolor='white', edgecolor='none', zorder=1)

    ax.set_title(title, fontsize=10, pad=25)
    #ax.text(0.5, -0.1, f"Global Total: {glob}", ha="center", transform=ax.transAxes)
    props = dict(boxstyle="round", facecolor='lightgray', alpha=0.5)
    text_box = ax.text(0.5, 1.07, f"Global Total: {glob}", ha="center", va="center", transform=ax.transAxes, bbox=props, fontsize=10)

    p = ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=0, vmax=masx)
    cbar = fig.colorbar(p, ax=ax, orientation=cborientation, fraction=fraction, pad=pad)
    cbar.set_label(clabel, labelpad=labelpad, fontsize=fontsize)

    return ax                 

################################################################

def define_diff_subplot(ax, data, lons, lats, cmap, cborientation, fraction, pad, labelpad, fontsize, title, glob, clabel):
    ax.coastlines(color='black')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.LAND, edgecolor='gray')

    ax.set_title(title)

    data_min, data_max = data.min(), data.max()

    if data_min == data_max:
        norm = mcolors.Normalize(vmin=data_min - 1, vmax=data_max + 1)
    else:
        abs_max = max(abs(0.25 * data_min), abs(0.25 * data_max))
        norm = mcolors.Normalize(vmin=-abs_max, vmax=abs_max)

    ax.set_title(title, fontsize=10, pad=22)
    #ax.text(0.5, -0.1, f"Global Total: {glob}", ha="center", transform=ax.transAxes)
    props = dict(boxstyle="round", facecolor='lightgray', alpha=0.5)
    text_box = ax.text(0.5, 1.07, f"Global Total: {glob}", ha="center", va="center", transform=ax.transAxes, bbox=props, fontsize=10)

    p = ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(), cmap=cmap, norm=norm)
    cbar = fig.colorbar(p, ax=ax, orientation=cborientation, fraction=fraction, pad=pad)
    cbar.set_label(clabel, labelpad=labelpad, fontsize=fontsize)
    return ax

################################################################

d = '/discover/nobackup/bkanzer/' 
d2 = '/discover/nobackup/bkanzer/nk_CCycle_E6obioF40'
d3 = '/discover/nobackup/projects/giss/prod_input_files/emis/BBURN_ALT/20240517/BBURN_GFED_4s/monthly'
d31 = '/discover/nobackup/projects/giss/prod_input_files/emis/BBURN_ALT/20240517/BBURN_GFED_4s/monthly/NAT'
d4 = '/discover/nobackup/projects/giss_ana/users/nkiang/aDATA/Biomass/Spawn2020/Spawn_Ent'
d5 = '/discover/nobackup/nkiang/DATA/Vegcov/V2HX2_EntGVSD_v1.1.2_16G_Spawn2020Sa_biomass_agb_2010_ann_pure.nc'

cmap = 'jet' #cividis, Spectral, jet
iyear = 2010
fyear = 2010
func_group = ["grass", "shrub", "tree"]

################################################################

srang = ["009", "010"]
trang = ["002", "004", "006", "007", "008", "009", "010"] #["001", "002", "003", "004", "005", "006", "007", "008", "009", "010"]
grang = ["011", "012", "013", "014"]

emis_G = np.zeros((90, 144),dtype=float)
emis_S = np.zeros((90, 144),dtype=float)    
emis_T = np.zeros((90, 144),dtype=float)    

for year in range(iyear, fyear + 1):
    filename = f"ANN{year}.aijnk_CCycle_E6obioF40.nc"
    filepath = os.path.join(d, filename)
 
    if os.path.exists(filepath):
        with nc.Dataset(filepath) as f:
            lats = f.variables['lat'][:]
            lons = f.variables['lon'][:]
            axyp = f.variables['axyp'][:]
            missing_val = -1.E-30
            nan_mat = np.full((len(lats), len(lons)), np.nan)
            zero_mat = np.zeros((len(lats), len(lons)), dtype=float)
            for vtype in func_group:
                if vtype == "grass":
                    BA_G = f.variables['BA_grass'][:]
                    BA_G = np.where(BA_G <= 0., zero_mat, BA_G)
                    emis_G = calculate_emis(vtype, BA_G, f, missing_val, nan_mat, srang, grang, trang)
                    #emis_G = np.where(emis_G <= -1000, zero_mat, emis_G)
                elif vtype == "shrub":
                    vtype = "shrub"
                    BA_S = f.variables['BA_shrub'][:]
                    BA_S = np.where(BA_S <= 0., zero_mat, BA_S)
                    emis_S = calculate_emis(vtype, BA_S, f, missing_val, nan_mat, srang, grang, trang)
                    #emis_S = np.where(emis_S <= -5000, zero_mat, emis_S) 
                if vtype == "tree":
                    BA_T = f.variables['BA_tree'][:]
                    BA_T = np.where(BA_T <= 0., zero_mat, BA_T)
                    emis_T = calculate_emis(vtype, BA_T, f, missing_val, nan_mat, srang, grang, trang)
                    #emis_T = np.where(emis_T <= -5000, zero_mat, emis_T)

            emis_tot = emis_G + emis_S + emis_T   # Units kgCm^(-2)yr^(-1)
            tot_emis = (np.nansum(emis_tot))*60*60*24*365
            tot_emis = format(tot_emis, '.3e')
            emis_min = emis_tot[:].min()
            offline_max = emis_tot[:].max()
            offline99 = np.percentile(emis_tot, 99.9)
        
    filename2 = f"ANN{year}.taijnk_CCycle_E6obioF40.nc"
    filepath2 = os.path.join(d2, filename2)

    if os.path.exists(filepath2):
        cO2v_data = []
        with nc.Dataset(filepath2) as f2:
            lats2 = f2.variables['lat'][:]
            lons2 = f2.variables['lon'][:]

            cO2v = f2.variables['CO2n_pyrE_src'][:]
            vmin2 = cO2v[:].min()
            vmax2 = cO2v[:].max()

            cO2v *= 1.0E-12

            axyp = f2.variables['axyp'][:]
            cO2v *= axyp
            cO2v *= (60.*60.*24.*365.)

            cO2v = np.where(cO2v <= 0., zero_mat, cO2v)

            for value in cO2v:
                cO2v_data.append(value)

        tot_cO2 = np.nansum(cO2v_data)  
        tot_cO2 = format(tot_cO2, '3e')
        diff_sdata = cO2v_data - emis_tot
        tot_diff1 = np.nansum(diff_sdata)
        tot_diff1 = format(tot_diff1, '3e')
        vmin3 = diff_sdata[:].min()
        #modelE_max = cO2v_data[:].max()
        modelE99 = np.percentile(cO2v_data, 99.9)

    filename3 = f"{year}.nc"
    filepath3 = os.path.join(d31, filename3)
    if os.path.exists(filepath3):
        missing_val2 = 0.
        ann_sum = np.zeros((len(lats2), len(lons2)), dtype=float)
        with nc.Dataset(filepath3) as f3:
            for k in range(12):
                
                GFED_data = f3.variables['CO2n'][k, :, :]  
                GFED_CO2 = GFED_data.reshape(90, 144)

                GFED_CO2  *= axyp
                GFED_CO2  *= (60.*60.*24.*365) 
                
                GFED_CO2 = np.where(GFED_CO2 <= 0., zero_mat, GFED_CO2)

                ann_sum += GFED_CO2
            totGFED = np.nansum(ann_sum)  
            totGFED = format(totGFED, '3e')
            ann_mean = ann_sum / 12.
            diff_srdata = cO2v_data - ann_mean
            tot_diff2 = np.nansum(diff_srdata)
            tot_diff2 = format(tot_diff2, '3e')
            diff_rsdata = emis_tot - ann_mean
            tot_diff3 = np.nansum(diff_rsdata)
            tot_diff3 = format(tot_diff3, '3e')
            GFED_max = GFED_CO2[:].max()
            GFED99 = np.percentile(ann_sum, 98.0)


    filename4 = f"V2HX2_EntGVSD_v1.1.2_16G_ann_pure_Spawn_AGB_2010.nc"
    filepath4 = os.path.join(d4, filename4)
    if os.path.exists(filepath4):
        for vtype in func_group:
            f = Dataset(r'/discover/nobackup/bkanzer/ANN2003.aijnk_CCycle_E6obioF40.nc', 'r')
            if vtype == "grass":
                BA_G = f.variables['BA_grass'][:]
                BA_G = np.where(BA_G <= 0., zero_mat, BA_G)
                spawnem_G = calculate_spawnemis(vtype, BA_G, zero_mat, d, f, srang, trang, grang)
                #spawnem_G = np.where(spawnem_G <= -1000, zero_mat, spawnem_G)
                spawnem_G *= axyp
                spawnem_G *= 44./12.
            elif vtype == "shrub":
                BA_S = f.variables['BA_shrub'][:]
                BA_S = np.where(BA_S <= 0., zero_mat, BA_S)
                spawnem_S = calculate_spawnemis(vtype, BA_S, zero_mat, d, f, srang, trang, grang)
                #spawnem_S = np.where(spawnem_S <= -5000, zero_mat, spawnem_S)
                spawnem_S *= axyp
                spawnem_S *= 44./12.
            elif vtype == "tree":
                BA_T = f.variables['BA_tree'][:]
                BA_T = np.where(BA_T <= 0., zero_mat, BA_T)
                spawnem_T = calculate_spawnemis(vtype, BA_T, zero_mat, d, f, srang, trang, grang)
                #spawnem_T = np.where(spawnem_S <= -5000, zero_mat, spawnem_T)
                spawnem_T *= axyp
                spawnem_T *= 44./12.
        spawnem_tot = spawnem_G + spawnem_S + spawnem_T
        tot_spawnem = (np.nansum(spawnem_tot))
        tot_spawnem = format(tot_spawnem, '3e')
        mod_spawn = cO2v_data - spawnem_tot
        mod_spawn = np.where(np.isnan(mod_spawn), zero_mat, mod_spawn)
        tot_diff4 = np.nansum(mod_spawn)
        tot_diff4 = format(tot_diff4, '3e')
        off_spawn = emis_tot - spawnem_tot
        tot_diff5 = np.nansum(off_spawn)
        tot_diff5 = format(tot_diff5, '3e')
        spawnem_max = spawnem_tot[:].max()
        spawnem99 = np.percentile(spawnem_tot, 99.0)
################################################################

tmin = min(0, emis_min, vmin2, vmin3)
#tmax = max(500, emis_max, vmax2, vmax3)

################################################################
frac = 0.05
pad = 0.025
lp = 2
fs = 10

#fig, axs = plt.subplots(3, 3, figsize=(18, 12), subplot_kw={'projection': ccrs.PlateCarree()})
#
#plt.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.9, hspace=0.16, wspace=0.1)
#
#define_subplot(axs[0, 0], emis_tot, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Offline Emissions {year}", tot_emis, 'Emissions [$kgCyr^{-1}$]', offline99)
#define_subplot(axs[0, 1], cO2v_data, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Online Emissions {year}", tot_cO2, 'Emissions [$kgyr^{-1}$]', modelE99)
#define_diff_subplot(axs[1, 0], diff_sdata, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Online and Offline", tot_diff1, 'Emissions [$kgCyr^{-1}-kgyr^{-1}$]')
##define_subplot(axs[0, 1], ann_mean, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"GFED4s Emissions {year}", totGFED, 'Emissions [$kgyr^{-1}$]', GFED99)
#define_diff_subplot(axs[2, 0], diff_srdata, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Online and GFED4s", tot_diff2, 'Emissions [$kgCyr^{-1}$]')
#define_diff_subplot(axs[1, 1], diff_rsdata, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Offline and GFED4s", tot_diff3, 'Emissions [$kgCyr^{-1}$]')
#define_subplot(axs[0, 2], spawnem_tot, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"S2020 Emission Data {year}", tot_spawnem, 'Emissions', spawnem99) #I need to fix this one
#define_diff_subplot(axs[1, 2], mod_spawn, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Online and S2020", tot_diff4, 'Emissions')
#define_diff_subplot(axs[2, 2], off_spawn, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Offline and S2020", tot_diff5, 'Emissions')
#
#fig.suptitle(f"Online, Offline, GFED4s, Spawn et al. 2020 Emission Data and Difference Plots", fontsize=20, fontweight='bold')
##plt.subplots_adjust(top=10)
#axs[2, 1].axis('off')
##plt.tight_layout(pad=5.0, w_pad=4.2, h_pad=70.0)
#
#plt.savefig('/home/bkanzer/PosterImage/3X3Plot.png')

fig, axs = plt.subplots(2, 4, figsize=(18, 12), subplot_kw={'projection': ccrs.PlateCarree()})

plt.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.9, hspace=0.16, wspace=0.1)

define_subplot(axs[0, 0], emis_tot, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Offline Emissions {year}", tot_emis, 'Emissions [$kgCyr^{-1}$]', offline99)
define_subplot(axs[0, 1], cO2v_data, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Online Emissions {year}", tot_cO2, 'Emissions [$kgyr^{-1}$]', modelE99)
define_diff_subplot(axs[0, 2], diff_sdata, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Online and Offline", tot_diff1, 'Emissions [$kgCyr^{-1}-kgyr^{-1}$]')
#define_subplot(axs[0, 1], ann_mean, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"GFED4s Emissions {year}", totGFED, 'Emissions [$kgyr^{-1}$]', GFED99)
define_diff_subplot(axs[0, 3], diff_srdata, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Online and GFED4s", tot_diff2, 'Emissions [$kgCyr^{-1}$]')
define_diff_subplot(axs[1, 0], diff_rsdata, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Offline and GFED4s", tot_diff3, 'Emissions [$kgCyr^{-1}$]')
define_subplot(axs[1, 1], spawnem_tot, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"S2020 Emission Data {year}", tot_spawnem, 'Emissions', spawnem99) #I need to fix this one
define_diff_subplot(axs[1, 2], mod_spawn, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Online and S2020", tot_diff4, 'Emissions')
define_diff_subplot(axs[1, 3], off_spawn, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Offline and S2020", tot_diff5, 'Emissions')

fig.suptitle(f"Online, Offline, GFED4s, Spawn et al. 2020 Emission Data and Difference Plots", fontsize=20, fontweight='bold')
plt.show()
#plt.subplots_adjust(top=10)
#axs[2, 1].axis('off')
#plt.tight_layout(pad=5.0, w_pad=4.2, h_pad=70.0)

#plt.savefig('/home/bkanzer/PosterImage/3X3Plot.png')

fig, axs = plt.subplots(1, 3, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})


define_subplot(axs[0], emis_tot, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Offline Emissions {year}", tot_emis, 'Emissions [$kgyr^{-1}$]', offline99)
define_subplot(axs[1], cO2v_data, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Online Emissions {year}", tot_cO2, 'Emissions [$kgyr^{-1}$]', modelE99)
define_diff_subplot(axs[2], diff_sdata, lons, lats, 'bwr', "horizontal", frac, pad, lp, fs, f"Difference of Online and Offline", tot_diff1, 'Emissions [$kgCyr^{-1}-kgyr^{-1}$]')
fig.suptitle(f"Online Offline Emission Data and Difference Plots 2010", fontsize=20, fontweight='bold', y=0.62)
plt.tight_layout(h_pad=0)
#plt.show()
#plt.savefig('/home/bkanzer/PosterImage/OnlineOffline.pdf')

#fig, axs = plt.subplots(4, 1, figsize=(50, 20), subplot_kw={'projection': ccrs.PlateCarree()})

#define_subplot(axs[0], spawnem_G, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Grass Emissions", 'Emissions')
#define_subplot(axs[1], spawnem_S, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Shrub Emissions", 'Emissions')   
#define_subplot(axs[2], spawnem_T, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Tree Emissions", 'Emissions')   
#define_subplot(axs[3], spawnem_tot, lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"Total Emissions", 'Emissions')   

#plt.show()
################################################################
d = Dataset(r'/discover/nobackup/bkanzer/ANN2003.aijnk_CCycle_E6obioF40.nc', 'r')
emis_file = Dataset('CC7.nc', 'w', format='NETCDF4')
for name, dimension in d.dimensions.items():
    emis_file.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

var_dim = d.variables['gpp'].dimensions
var_dtype = d.variables['gpp'].datatype

totemis_var = emis_file.createVariable('emis_tot', var_dtype, var_dim)
emisG_var = emis_file.createVariable('emis_G', var_dtype, var_dim)
emisS_var = emis_file.createVariable('emis_S', var_dtype, var_dim)
emisT_var = emis_file.createVariable('emis_T', var_dtype, var_dim)

totemis_var[:] = emis_tot
emisG_var[:] = emis_G
emisS_var[:] = emis_S
emisT_var[:] = emis_T

d.close
