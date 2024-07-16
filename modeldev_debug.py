import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import sys

#devtype = 'setuphour'
devtype = 'month'
#devtype = 'pyrEsource'
if devtype == 'month':
    #pyrEsource 
    month='JAN1950'
    outputpath = '/discover/nobackup/kmezuman/plots/fire_repository/Development/setuphour_bug/'
    simlist = ['E6TomaF40pyrEold25', 'E6TomaF40pyrE25']
    index_pairs = [(0, 1)]
elif devtype == 'pyrEsource':
    outputpath = '/discover/nobackup/kmezuman/plots/fire_repository/Development/pyrEsource/'
    month='ANN1996'
    simlist = ['E6TomaF40intpyrEtest2', 'E6TomaF40intpyrEtest']
    index_pairs = [(0, 1)]
elif devtype == 'setuphour':
    month='PARTIAL'
    outputpath = '/discover/nobackup/kmezuman/plots/fire_repository/Development/setuphour_bug/Partial/'
    simlist = ['E6TomaF40pyrEold241', 'E6TomaF40pyrEold25', 'E6TomaF40pyrE241', 'E6TomaF40pyrE25']
    index_pairs = [(0, 1), (2, 3), (2, 0), (3, 1)]
else:
    sys.exit(f'incompatibel type: {devtype}')

simpath = '/discover/nobackup/kmezuman/'
var_list = ['fireCount']#, 'BA_tree', 'BA_grass', 'BA_shrub', 'FLAMM', 'f_ignCG', 'f_ignHUMAN', 'CO_biomass_src']
hemis_glob_ind = 2
data = {}

# Loop over relevant variables
for s, sim in enumerate(simlist):
    for v, var in enumerate(var_list):
    # Set ftype depending on the var 
        if var.endswith('biomass_src') or var.endswith('pyrE_src'):
            ftype = 'taij'
        else:
            ftype = 'aij'

        filename = f"{simpath}{sim}/{month}.{ftype}{sim}.nc"
        try:
            dataset = nc.Dataset(filename, 'r')

            if v == 0 and s == 0:
                lat = dataset.variables['lat'][:]
                lon = dataset.variables['lon'][:]
                hemis_nan_array = np.full(3, np.nan)
                nan_mat = np.full((len(lat), len(lon)), np.nan)
            # Read hemis values
            hemis_var = dataset.variables[var + '_hemis']
            # Read 2D values
            two_d_var = dataset.variables[var]
            # Handle missing values for hemis_var 
            fill_value = hemis_var._FillValue if '_FillValue' in hemis_var.ncattrs() else None 
            missing_value = hemis_var.missing_value if 'missing_value' in hemis_var.ncattrs() else None 
            hemis_val = hemis_var[:] 
            if missing_value is not None: 
                hemis_val = np.where(hemis_val <= missing_value, hemis_nan_array, hemis_val) 
            # Handle missing values for two_d_var 
            missing_value_2d = two_d_var.missing_value if 'missing_value' in two_d_var.ncattrs() else None 
            two_d_val = two_d_var[:] 
            if missing_value_2d is not None: 
                two_d_val = np.where(two_d_val <= missing_value_2d, nan_mat, two_d_val) 
            data[(var, s)] = {"hemis": hemis_val[hemis_glob_ind], "2D": two_d_val}
        except FileNotFoundError:
            print(f"File {filename} not found.")
        except Exception as error:
            print(f"An exception occurred while reading from {filename}: {str(var)}. {error}")
for key, value in data.items():
    var_name, sim_index = key
    hemis_shape = value["hemis"].shape
    two_d_shape = value["2D"].shape
    #print(f"Variable: {var_name}, Simuation Index: {sim_index}")
    #print(f"  Hemis Shape: {hemis_shape}")
    #print(f"  2D Shape: {two_d_shape}")

# Calculate differences
hemis_diffs_dict = {var: [] for var in var_list}
simnames = []
for i, j in index_pairs:
    if i-j > 1:
        simname = f"({simlist[i]} - {simlist[j]})/{simlist[j]}"
    else:
        simname = f"{simlist[i]} - {simlist[j]}"
    print(simname)
    simnames.append(simname)
    for v, var in enumerate(var_list):
        if i-j > 1:
            hemis_diff = (data[(var, i)]["hemis"] - data[(var, j)]["hemis"])/data[(var, j)]["hemis"]
        else:
            hemis_diff = data[(var, i)]["hemis"] - data[(var, j)]["hemis"]
        hemis_diffs_dict[var].append(hemis_diff)

# Print and save the difference into a text file
df = pd.DataFrame(hemis_diffs_dict, index = simnames)
formatted_table = df.to_string()
output_file_path = (f"{outputpath}fire_total_diff.txt")
with open(output_file_path, 'w') as file:
    file.write(formatted_table)
print(f"Differences saved in {output_file_path}")

#Plot 2D differences
for var in var_list:
    #plt.subplots_adjust(hspace=0.001)
    if devtype == 'month':
        plots = [
            (0, 0, (var, 0), simlist[0]), # A 
            (0, 1, (var, 1), simlist[1]), # B 
            (0, 2, (0, 1), f"{simlist[0]} - {simlist[1]}"), # A-B 
            ]
        fig, axes = plt.subplots(nrows=1, ncols=len(plots), figsize=(20, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    elif devtype == 'setuphour':
        plots = [
            (0, 0, (var, 0), simlist[0]), # A 
            (0, 1, (var, 1), simlist[1]), # B 
            (0, 2, (0, 1), f"{simlist[0]} - {simlist[1]}"), # A-B 
            (1, 0, (var, 2), simlist[2]), # A' 
            (1, 1, (var, 3), simlist[3]), # B' 
            (1, 2, (2, 3), f"{simlist[2]} - {simlist[3]}"), # A'-B' 
            (2, 0, (0, 2), f"{simlist[0]} - {simlist[2]}"), # A-A' 
            (2, 1, (1, 3), f"{simlist[1]} - {simlist[3]}"), # B-B'
             ] 
        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(20, 15), subplot_kw={'projection': ccrs.PlateCarree()})
        fig.delaxes(axes[2,2])
    else:
        sys.exit(f'incompatibel type: {devtype}')
    axes=axes.flatten()[:len(plots)]
    for i, (row, col, plot_info, label) in enumerate(plots): 
        ax = axes[i] 
        ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

        if isinstance(plot_info[0], str): # Plotting actual values A, B, A', B' 
            var, idx = plot_info
            plot_data = data[(var, idx)]["2D"] 
            cmap = 'jet'
            vmax = 0.7 * np.nanmax(plot_data)
            vmin = 0.
        else: # Plotting differences 
            i, j = plot_info          
            plot_data = data[(var, i)]["2D"] - data[(var, j)]["2D"] 
            cmap='bwr'
            max_diff = 0.7 * np.nanmax(np.abs(plot_data))
            vmax = max_diff
            vmin = -max_diff

        total = np.nansum(plot_data)
        #total_str = f"{total:.2f}"
        ax.set_title(f"{var}: {label}\nTotal: {total}", y=1.10)
        longitudes, latitudes = np.meshgrid(lon, lat)
        plot = ax.pcolormesh(longitudes, latitudes, plot_data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(plot, ax=ax, orientation='horizontal',fraction=0.046, pad=0.08)
        cbar.set_label('units',labelpad=5,fontsize=10)
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS)
        #ax.add_feature(cfeature.LAND)
        #ax.add_feature(cfeature.OCEAN)
        gl = ax.gridlines(Visible=None,draw_labels=True)
        #gl.draw_gridlines(False)
        #gl.top_labele = False
        #gl.right_labels = False
        #gl.xlabel_style = {'size': 10, 'color': 'gray'} 
        #gl.ylabel_style = {'size': 10, 'color': 'gray'}
        #fig.colorbar(plot, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)

    fig.savefig(f"{outputpath}{var}_eval.eps")
