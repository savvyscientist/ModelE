import os
import numpy as np
import argparse
import xarray as xr


def read_nc_list(path):

    for root, dirs, files in os.walk(path):

        sorted_files = sorted(files)
        nc_files_list = [os.path.join(root, file) for file in sorted_files if file.split(".")[1] == "taijlh1E33oma_ai"]

    return nc_files_list


def main(args):

    ds = xr.open_mfdataset(read_nc_list(args.root_path))
    
    print("Available dimensions:")
    for key, value in ds.sizes.items():
        print("Dimension name:", key, "\t Dimeniton length:", value)
    
    weights = np.cos(np.deg2rad(ds.lat))
    weights.name = "weights"
    
    # Ask the user to select one or more dimensions, separated by commas
    selected_dims = input("Please enter the dimensions you want to take the average over (comma-separated): ")

    # Convert the user input into a list of dimensions
    selected_dims_list = [dim.strip() for dim in selected_dims.split(',')]
    
    # Validate the selected dimensions
    dims = list(ds.sizes.keys())
    invalid_dims = [dim for dim in selected_dims_list if dim not in dims]

    if invalid_dims:
        print(f"Invalid dimension(s) selected: {invalid_dims}. Please try again.")

    else:
        # Take the average over the selected dimensions
        averaged_ds = ds[args.species].mean(dim=selected_dims_list)
        print(averaged_ds.shape)
         
        
def get_arguments():
    parser = argparse.ArgumentParser(description="Analyzing diagnostics of ModelE in NetCDF format.")

    parser.add_argument("--root-path", type=str, default="/home/serfani/serfani_data1/E33OMA", 
        help="Root path for NetCDF file.")
    parser.add_argument("--species", type=str, default="BCB", 
        help="Aerosol species.")
    

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_arguments()
    main(args)

