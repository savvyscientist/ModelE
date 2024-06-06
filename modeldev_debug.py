import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd

var_list = ['fireCount', 'BA_tree', 'BA_grass', 'BA_shrub', 'FLAMM', 'f_ignCG', 'f_ignHUMAN', 'CO_biomass_src']
simpath = '/discover/nobackup/kmezuman/'
outputpath = '/discover/nobackup/kmezuman/plots/fire_repository/Development/setuphour_bug/'
simlist = ['E6TomaF40pyrEold241', 'E6TomaF40pyrEold25', 'E6TomaF40pyrE241', 'E6TomaF40pyrE25']
hemis_glob_ind = 2
data = {}

# Loop over relevant variables
for s, sim in enumerate(simlist):
    for v, var in enumerate(var_list):
    # Set ftype depending on the var 
        if var.endswith('biomass_src'):
            ftype = 'taij'
        else:
            ftype = 'aij'

        filename = f"{simpath}{sim}/PARTIAL.{ftype}{sim}.nc"

        try:
            dataset = nc.Dataset(filename, 'r')

            if v == 0 and s == 0:
                lat = dataset.variables['lat'][:]
                lon = dataset.variables['lon'][:]
            # Read hemis values
            hemis_var = dataset.variables[var + '_hemis']
            # Read 2D values
            two_d_var = dataset.variables[var]
            # Handle missing values 
            fill_value = hemis_var._FillValue if '_FillValue' in hemis_var.ncattrs() else np.nan 
            hemis_val = np.where(hemis_var[:] == fill_value, np.nan, hemis_var[:]) 
            fill_value_2d = two_d_var._FillValue if '_FillValue' in two_d_var.ncattrs() else np.nan 
            two_d_val = np.where(two_d_var[:] == fill_value_2d, np.nan, two_d_var[:]) 
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
index_pairs = [(0, 1), (2, 3), (0, 2), (1, 3)]
hemis_diffs_dict = {var: [] for var in var_list}
simnames = []
for i, j in index_pairs:
    simname = f"{simlist[i]} - {simlist[j]}"
    simnames.append(simname)
    for v, var in enumerate(var_list):
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
    fig, axes = plt.subplots(nrows=len(index_pairs), ncols=1, figsize=(20, 30), subplot_kw={'projection': ccrs.PlateCarree()})
    if len(index_pairs) == 1:
        axes = [axes]
    ind = 0
    for i, j in index_pairs:
        two_d_diff = data[(var, i)]["2D"] - data[(var, j)]["2D"]
        ax = axes[ind]
        ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
        ax.set_title(f"{var}: {simname}", y=1.05)
        longitudes, latitudes = np.meshgrid(lon, lat)
        plot = ax.pcolormesh(longitudes, latitudes, two_d_diff, transform=ccrs.PlateCarree(), cmap='viridis')
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.OCEAN)
        gl = ax.gridlines(draw_labels=True)
        gl.top_labele = False
        gl.right_labels = False
        fig.colorbar(plot, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
        ind += 1

    plt.savefig(f"{outputpath}{var}_diff_plots.eps")
    plt.close(fig)
