import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import sys
import re

#devtype = 'setuphour'
#devtype = 'month'
devtype = 'pyrEsource'
if devtype == 'month':
    month = 'JAN1950'
    outputpath = '/discover/nobackup/kmezuman/plots/fire_repository/Development/setuphour_bug/'
    simlist = ['E6TomaF40pyrEold25', 'E6TomaF40pyrE25']
    index_pairs = [(0, 1)]
elif devtype == 'pyrEsource':
    outputpath = '/discover/nobackup/kmezuman/plots/fire_repository/Development/pyrEsource/'
    month = 'ANN1996'
    simlist = ['E6TomaF40intpyrEtest2', 'E6TomaF40intpyrEtest']
    index_pairs = [(0, 1)]
elif devtype == 'setuphour':
    month = 'PARTIAL'
    outputpath = '/discover/nobackup/kmezuman/plots/fire_repository/Development/setuphour_bug/Partial/'
    simlist = ['E6TomaF40pyrEold241', 'E6TomaF40pyrEold25', 'E6TomaF40pyrE241', 'E6TomaF40pyrE25']
    index_pairs = [(0, 1), (2, 3), (2, 0), (3, 1)]
else:
    sys.exit(f'incompatibel type: {devtype}')

simpath = '/discover/nobackup/kmezuman/'

# Define the baseline variables
#base_variables = ["FLASH", "CtoG", "FLAMM", "fireCount", "f_ignCG", "f_ignHUMAN", "BA_tree", "BA_shrub", "BA_grass"]
base_variables = ["fireCount","FLAMM"]
s_to_yr = 60.*60.*24.*365.
kg_to_Tg = 1.E-9

# Define the tracer sources and tracers
tracer_sources = ["AGRI_CMIP6_src", "PEAT_CMIP6_src", "biomass_src", "pyrE_src"]
#tracers = ["CO", "BCB"]
tracers = ["CO", "BCB", "OCB", "Alkenes", "Paraffin", "SO2", "NOx", "NH3"]
include_tracers = True 

# Define the carbon cycle tracer
ccycle = "CO2n"
include_ccycle = True

intch4 = 'CH4'
include_intch4 = False

# Initialize the variable list with baseline variables
var_list = base_variables.copy()

# Add tracer-related variables
if include_tracers:
    for tracer in tracers:
        for source in tracer_sources:
            var_list.append(f"{tracer}_{source}")

# Add carbon cycle-related variables
if include_ccycle:
    for source in tracer_sources:
        var_list.append(f"{ccycle}_{source}")

# Add methane related variables
if include_intch4:
    for source in tracer_sources:
        var_list.append(f"{intch4}_{source}")

hemis_glob_ind = 2
data = {}

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

# A function that changes a key's name in a dictionary
def change_dict_key(d, old_key, new_key, default_value=None):
    d[new_key] = d.pop(old_key, default_value)

# Loop over simulations and variables to process data
for s, sim in enumerate(simlist):
    for v, var in enumerate(var_list):
        # Determine the file type based on variable name
        if var.endswith('_src'):
            ftype = 'taij'
        else:
            ftype = 'aij'

        filename = f"{simpath}{sim}/{month}.{ftype}{sim}.nc"
        try:
            dataset = nc.Dataset(filename, 'r')

            if v == 0 and s == 0:
                lat = dataset.variables['lat'][:]
                lon = dataset.variables['lon'][:]
                grid_cell_area = dataset.variables['axyp'][:]
                #hemis_nan_array = np.full(1, np.nan)
                hemis_nan_array = np.full(1, 0.)
                #nan_mat = np.full((len(lat), len(lon)), np.nan)
                nan_mat = np.full((len(lat), len(lon)), 0.)

            # Extract units of the variable and the scaling factor
            units = dataset.variables[var].units
            scaling_factor, unit = extract_scaling_factor(units)
            print(var,scaling_factor,units)

            #Calculate global mean values
            hemis_var = dataset.variables[var + '_hemis']
            hemis_var = hemis_var[hemis_glob_ind]
            # Handle missing values for hemis_var 
            #the line commented below has a bug
            #missing_value = hemis_var.missing_value #if 'missing_value' in hemis_var.ncattrs() else None 
            missing_value = -1.E+30
            if missing_value is not None: 
                hemis_var = np.where(hemis_var <= missing_value, hemis_nan_array, hemis_var) 
            hemis_var *= scaling_factor

            # Read 2D values
            two_d_var = dataset.variables[var]
            # Handle missing values for two_d_var 
            #the line commented below has a bug
            #missing_value_2d = two_d_var.missing_value if 'missing_value' in two_d_var.ncattrs() else None 
            missing_value_2d = missing_value
            two_d_val = two_d_var[:] 
            if missing_value_2d is not None: 
                two_d_val = np.where(two_d_val <= missing_value_2d, nan_mat, two_d_val) 
            two_d_val *= scaling_factor

            #Calculate global total values
            if "m^-2" in unit:
                 total = np.nansum(two_d_val * grid_cell_area)
                 total_unit = unit.replace("m^-2","")
            elif "m-2" in unit:
                 total = np.nansum(two_d_val * grid_cell_area)
                 total_unit = re.sub(r'\bm\-2\b', '', unit).strip()
            else:
                 total = np.nansum(two_d_val)
                 total_unit = unit

            if "s-1" in unit:
                 total *= s_to_yr
                 total_unit = total_unit.replace("s-1","yr-1").strip()

            if "kg" in unit:
                 total *= kg_to_Tg 
                 total_unit = total_unit.replace("kg","Tg").strip()


            # Store data in dictionary
            data[(var, s)] = {
                    "mean": hemis_var,
                    "total": total, 
                    "total_units": total_unit,
                    "2D": two_d_val,
                    "units": unit
            }

            #print(hemis_var.shape, two_d_val.shape)
            dataset.close()

        except FileNotFoundError:
            print(f"File not found: {filename}")

        except Exception as error:
            print(f"An exception occurred while reading from {filename}: {str(var)}. {error}")
            data[(var, s)] = {"hemis": np.zeros(1,dtype=float), "2D": np.zeros((len(lat), len(lon)))}
            # Store data in dictionary
            data[(var, s)] = {
                    "mean": np.zeros(1,dtype=float),
                    "total": np.zeros(1,dtype=float), 
                    "total_units": [],
                    "2D": np.zeros((len(lat), len(lon))),
                    "units": []
            }

#for key, value in data.items():
#    var_name, sim_index = key
#    two_d_shape = value["2D"].shape
#    print(f"Variable: {var_name}, Simuation Index: {sim_index}")
#    print(f"  2D Shape: {two_d_shape}")

# Calculate differences
global_diffs_dict = {var: [] for var in var_list}
simnames = []
for i, j in index_pairs:
    if i-j > 1:
        simname = f"({simlist[i]} - {simlist[j]})/{simlist[j]}"
    else:
        simname = f"{simlist[i]} - {simlist[j]}"
    #print(simname)
    simnames.append(simname)
    for v, var in enumerate(var_list):
        if i-j > 1:
            global_diff = (data[(var, i)]["total"] - data[(var, j)]["total"])/data[(var, j)]["total"]
        else:
            global_diff = data[(var, i)]["total"] - data[(var, j)]["total"]
        global_diffs_dict[var].append(global_diff)

#print(data.items())
#print(data.keys())
#print(data.values())
# Print and save the totals into a text file
totals_dict = {}
# Extract the "total" values for each variable and simulation 
for (var, sim_index), values in data.items():
    #print(var,sim_index,simlist[sim_index],values["total"])
    var_with_unit = f"{var} [{values['total_units']}]"
    #print(var_with_unit)
    if var not in totals_dict:
        totals_dict[var] = {}
    totals_dict[var][simlist[sim_index]] = values["total"]
#the construction of totals_dict has a bug
#add a part that checks if 'total_units' are [] and if so rename the dictionary
#key:
#    if values['total_units']
#    change_dict_key(totals_dict, var, var_with_unit, default_value=None):

# Create a DataFrame from the totals_dict 
df_totals = pd.DataFrame(totals_dict) 
# Transpose the DF to have simulations as a columns and variables as rows
df_totals = df_totals.T
# Print the DataFrame to check the result 
print(df_totals)
formatted_table = df_totals.to_string()
output_file_path = (f"{outputpath}fire_total.txt")
with open(output_file_path, 'w') as file:
    file.write(formatted_table)
print(f"Totals saved in {output_file_path}")

# Print and save the difference into a text file
df = pd.DataFrame(global_diffs_dict, index = simnames)
formatted_table = df.to_string()
output_file_path = (f"{outputpath}fire_total_diff.txt")
with open(output_file_path, 'w') as file:
    file.write(formatted_table)
print(f"Differences saved in {output_file_path}")

#Plot 2D differences
for var in var_list:
    #plt.subplots_adjust(hspace=0.001)
    if devtype == 'month' or devtype == 'pyrEsource':
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
    axes = axes.flatten()[:len(plots)]
    for i, (row, col, plot_info, label) in enumerate(plots): 
        ax = axes[i] 
        ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())

        #calculate the 95th percentile of plot_data and use that as the colorbar max

        if isinstance(plot_info[0], str): # Plotting actual values A, B, A', B' 
            var, idx = plot_info
            plot_data = data[(var, idx)]["2D"] 
            cmap = 'jet'
            vmax = np.percentile(plot_data, 97)  # Use the 95th percentile for colorbar max
            vmin = 0.
        else: # Plotting differences 
            i, j = plot_info          
            plot_data = data[(var, i)]["2D"] - data[(var, j)]["2D"] 
            cmap='bwr'
            max_diff = np.percentile(np.abs(plot_data), 97)  # Use the 95th percentile for colorbar max
            vmax = max_diff
            vmin = -max_diff

        ax.set_title(f"{var}: {label}\nTotal: {data[(var, idx)]['total']} {data[(var, idx)]['total_units']}", y=1.10)
        longitudes, latitudes = np.meshgrid(lon, lat)
        plot = ax.pcolormesh(longitudes, latitudes, plot_data, transform=ccrs.PlateCarree(), cmap=cmap, vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(plot, ax=ax, orientation='horizontal',fraction=0.046, pad=0.08)
        #fix units
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
    plt.close(fig)
