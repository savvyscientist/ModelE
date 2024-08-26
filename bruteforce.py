import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import subprocess 
import numpy as np
import sys
import netCDF4 as nc
'''
I believe this code calculates, and plots, the vegetation fraction for each PFT
'''
def define_subplot(ax, data, lons, lats, cmap, cborientation, fraction, pad, labelpad, fontsize, title, clabel):
    ax.coastlines(color='black')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.LAND, edgecolor='gray')

    ax.set_title(title)
    p = ax.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(), cmap=cmap)
    #cbar = fig.colorbar(p, ax=ax, orientation=cborientation, fraction=fraction, pad=pad)
    #cbar.set_label(" ", labelpad=labelpad, fontsize=fontsize)

    return ax


iyear = 2000
fyear = 2015

frac = 0.05
pad = 0.001
lp = 5
fs = 1

d = '/discover/nobackup/bkanzer'
d2 = '/discover/nobackup/bkanzer/ANNVar_Figures'

for year in range(iyear, fyear + 1):
    filename = f"ANN{year}.aijnk_CCycle_E6obioF40.nc"
    filepath = os.path.join(d, filename)

    if os.path.exists(filepath):
        with nc.Dataset(filepath) as f:
            lats = f.variables['lat'][:]
            lons = f.variables['lon'][:]
            nan_mat = np.full((len(lats), len(lons)), np.nan)
            missingval = 0.
            #for g in range(1, 17):
                #print(f"ra0010{str(g).zfill(2)}")
            #sys.exit()
            ra_vars = [f.variables[f"ra0010{str(i).zfill(2)}"][:] for i in range(1, 17)]
            #ra_vars = [np.where(ra <= missingval, nan_mat, ra) for ra in ra_vars]

            fig, axs = plt.subplots(4, 4, figsize=(30, 20), subplot_kw={'projection': ccrs.PlateCarree()})
            lon, lat = np.meshgrid(lons, lats)
            cmap = 'jet'
            frac, pad, lp, fs = 0.036, 0.06, 0.1, 10

            for i, ax in enumerate(axs.flat):
                #print(i)
                if i == 14 or i == 15:
                    pass
                else:
                    define_subplot(ax, ra_vars[i], lons, lats, cmap, "horizontal", frac, pad, lp, fs, f"ra0010{str(i + 1).zfill(2)}", f"{i + 1}")

            #plt.show()
        filesave = os.path.join(d2, f"Var_Plot{year}.eps")
        plt.savefig(filesave)
