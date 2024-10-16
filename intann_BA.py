import h5py
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def read_gfed4s(file_paths):
    """
    Reads multiple HDF5 files using h5py, calculates the annual burned area,
    and returns the data as xarray.DataArray.
    """
    burned_fraction_list = []

    for file_path in file_paths:
        # Open the HDF5 file using h5py
        with h5py.File(file_path, 'r') as h5file:
            # Access grid_cell_area using the method suggested
            grid_cell_area = h5file['ancill']['grid_cell_area'][:]
            
            # Load lat and lon for constructing the xarray dataset
            lat = h5file['lat'][:]
            lon = h5file['lon'][:]
            
            # Sum burned fraction over all months
            annual_burned_fraction = np.zeros_like(grid_cell_area)
            for month in range(1, 13):
                month_burned_fraction = h5file[f'burned_area/{month:02d}/burned_fraction'][:]
                annual_burned_fraction += month_burned_fraction

            # Calculate total burned area
            total_burned_area = annual_burned_fraction * grid_cell_area
            burned_fraction_list.append(total_burned_area)

    # Convert the list to xarray.DataArray for further processing
    total_burned_area_all_years = xr.DataArray(
        burned_fraction_list, 
        dims=['year', 'lat', 'lon'], 
        coords={'lat': lat, 'lon': lon}
    )

    return total_burned_area_all_years

def read_ModelEBA(startyear, endyear, simname, ModelE_path):
    """
    Reads ModelE BA data (BA_tree, BA_shrub, BA_grass) for the given year range, sums them to calculate 
    modelE_BA, and returns the annual sum for each year.

    Parameters:
    startyear (int): The starting year.
    endyear (int): The ending year.
    simname (str): The simulation name to match the file format (default 'nk_CCycle_E6obioF40').
    ModelE_path (str): The directory containing ModelE output.

    Returns:
    np.ndarray: A 2D array (year, modelE_BA), where modelE_BA is the sum of BA_tree, BA_shrub, and BA_grass.
    """

    # Generate a list of file paths for the given year range
    file_paths = [f'{ModelE_path}ANN{year}.aij{simname}.nc' for year in range(startyear, endyear + 1)]
    
    # Use xarray's open_mfdataset to open all files and concatenate along the 'time' dimension (which represents year)
    ds = xr.open_mfdataset(file_paths, combine='by_coords')
    
    # Sum BA_tree, BA_shrub, and BA_grass to get modelE_BA
    modelE_BA = ds['BA_tree'] + ds['BA_shrub'] + ds['BA_grass']
    
    return modelE_BA

def intann_BA_xarray(startyear, endyear, GFED_path, ModelE_path, simname):
    """
    Calculates the decade mean burned area (BA) and the interannual variability of BA
    from 2002 to 2012 using read_gfed4s and read_ModelEBA.

    Parameters:
    startyear (int): The start year of the period (default 2002).
    endyear (int): The end year of the period (default 2012).
    GFED_path (str): The directory containing the GFED4s files.
    ModelE_path (str): The directory containing ModelE output.
    simname (str): The simulation name for ModelE data.

    Returns:
    tuple: A tuple containing:
      - decade_mean_ba (xarray.DataArray): The mean burned area over the decade (lat, lon array).
      - ba_per_year (np.ndarray): A 2D array with columns (year, totalBA), where totalBA is the sum of burned area for that year.
      - modelE_BA_per_year (np.ndarray): A 2D array with columns (year, modelE_BA), where modelE_BA is the total burned area from ModelE.
    """

    # Call read_gfed4s to load GFED4s data
    file_paths = [f'{GFED_path}GFED4.1s_{year}.hdf5' for year in range(startyear, endyear + 1)]
    total_burned_area_all_years = read_gfed4s(file_paths)

    # Calculate the mean burned area over the decade
    decade_mean_ba = total_burned_area_all_years.mean(dim='year')

    # Calculate total burned area for each year from GFED4s data
    total_ba_per_year = total_burned_area_all_years.sum(dim=['phony_dim_0', 'phony_dim_1']).values
    years = np.arange(startyear, endyear + 1)
    ba_per_year = np.column_stack((years, total_ba_per_year))

    # Call read_ModelEBA to load and process ModelE data
    modelE_BA_all_years = read_ModelEBA(startyear, endyear, simname, ModelE_path)

    # Calculate the mean burned area over the decade (ModelE)
    decade_mean_modelEba = modelE_BA_all_years.mean(dim='time')

    # Calculate total burned area for each year from ModelE data
    total_modelE_BA_per_year = modelE_BA_all_years.sum(dim=['lat', 'lon']).values
    modelE_BA_per_year = np.column_stack((years, total_modelE_BA_per_year))

    return decade_mean_ba, decade_mean_modelEba, ba_per_year, modelE_BA_per_year 

def map_plot(decade_mean_ba):
    """
    Plots the decadal mean burned area of both GFED and ModelE side by side.
    
    Parameters:
    decade_mean_ba (xarray.DataArray): The decadal mean burned area (lat, lon array).
    decade_mean_modelEba (xarray.DataArray): The decadal mean burned area from ModelE(lat, lon array).
    """
   
    # Plot side by side maps for GFED and ModelE
    fig, axes = plt.subplots(ncols=2, figsize=(14, 6))

    # GFED decadal mean map
    decade_mean_ba.plot(ax=axes[0], cmap='YlOrRd')
    axes[0].set_title('GFED Decadal Mean Burned Area (2002-2012)')

    # ModelE decadal mean map 
    decade_mean_modelEba.plot(ax=axes[1], cmap='YlOrRd')
    axes[1].set_title('ModelE Decadal Mean Burned Area (2002-2012)')

    plt.show()

def time_series_plot(ba_per_year, modelE_BA_per_year):
    """
    Plots the total burned area as a function of year for both GFED and ModelE data.
    
    Parameters:
    ba_per_year (np.ndarray): A 2D array with columns (year, totalBA), where totalBA is the sum of burned area for that year.
    modelE_BA_per_year (np.ndarray): A 2D array with columns (year, modelE_BA), where modelE_BA is the sum of burned area for that year.
    """
    
    # Extract years and total burned area for both GFED and ModelE
    years_gfed = ba_per_year[:, 0]
    total_ba_gfed = ba_per_year[:, 1]
    
    years_modelE = modelE_BA_per_year[:, 0]
    total_ba_modelE = modelE_BA_per_year[:, 1]
    
    # Plot the time series of total burned area for both GFED and ModelE
    plt.figure(figsize=(10, 6))
    plt.plot(years_gfed, total_ba_gfed, marker='o', linestyle='-', color='b', label='GFED BA')
    plt.plot(years_modelE, total_ba_modelE, marker='x', linestyle='-', color='r', label='ModelE BA')
    plt.title('Total Burned Area Over Time (2002-2012)')
    plt.xlabel('Year')
    plt.ylabel('Total Burned Area (BA)')
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    # Example usage with test parameters
    startyear = 2002
    endyear = 2012
    GFED_path = '/discover/nobackup/projects/giss_ana/users/kmezuman/gfed4s/'  
    ModelE_path = '/discover/nobackup/kmezuman/nk_CCycle_E6obioF40/'
    simname = 'nk_CCycle_E6obioF40'

    # Call intann_BA_xarray to calculate decadal mean BA and interannual variability
    decade_mean_ba, ba_per_year, modelE_BA_per_year = intann_BA_xarray(startyear, endyear, GFED_path, ModelE_path, simname)

    # Plot the decadal mean burned area
    map_plot(decade_mean_ba, decade_mean_modelEba)

    # Plot the time series of burned area for GFED and ModelE
    time_series_plot(ba_per_year, modelE_BA_per_year)

if __name__ == '__main__':
    main()
