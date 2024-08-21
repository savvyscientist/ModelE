import xarray as xr
import numpy as np

def conservative_regrid_manual(burned_area_data, target_lat_points, target_lon_points):
    """
    Manually perform conservative regridding from a high-resolution grid to a lower-resolution grid.
    
    Parameters:
    - burned_area_data: xarray.DataArray
        The high-resolution burned area data to be regridded.
    - target_lat_points: int
        The number of latitude points in the target resolution.
    - target_lon_points: int
        The number of longitude points in the target resolution.
        
    Returns:
    - regridded_data: xarray.DataArray
        The regridded burned area data.
    """
    # Define original grid
    lat_old = burned_area_data.lat
    lon_old = burned_area_data.lon

    # Create target grid
    lat_new = np.linspace(lat_old.min(), lat_old.max(), target_lat_points)
    lon_new = np.linspace(lon_old.min(), lon_old.max(), target_lon_points)
    
    # Initialize the new regridded array
    regridded_data_values = np.zeros((target_lat_points, target_lon_points))
    
    # Calculate the grid cell size for both grids
    lat_old_cell_size = (lat_old[1] - lat_old[0]).item()
    lon_old_cell_size = (lon_old[1] - lon_old[0]).item()
    lat_new_cell_size = (lat_new[1] - lat_new[0])
    lon_new_cell_size = (lon_new[1] - lon_new[0])

    # Loop over each cell in the target grid
    for i in range(target_lat_points):
        for j in range(target_lon_points):
            # Find the bounds of the new grid cell
            lat_min = lat_new[i] - lat_new_cell_size / 2
            lat_max = lat_new[i] + lat_new_cell_size / 2
            lon_min = lon_new[j] - lon_new_cell_size / 2
            lon_max = lon_new[j] + lon_new_cell_size / 2

            # Find which cells in the original grid contribute to this new cell
            lat_indices = np.where((lat_old >= lat_min) & (lat_old < lat_max))[0]
            lon_indices = np.where((lon_old >= lon_min) & (lon_old < lon_max))[0]

            # Accumulate the values from the original grid into the new grid cell
            for lat_index in lat_indices:
                for lon_index in lon_indices:
                    overlap_area = lat_old_cell_size * lon_old_cell_size  # Approximate overlap area
                    regridded_data_values[i, j] += burned_area_data[lat_index, lon_index].item() * overlap_area

            # Normalize by the area of the new cell to approximate conservation
            regridded_data_values[i, j] /= lat_new_cell_size * lon_new_cell_size
    
    # Create the new xarray DataArray for the regridded data
    regridded_data = xr.DataArray(
        regridded_data_values,
        coords=[lat_new, lon_new],
        dims=["lat", "lon"]
    )
    
    return regridded_data

# Example usage
lat_old = np.linspace(-89.5, 89.5, 180)
lon_old = np.linspace(-179.5, 179.5, 360)
burned_area_data = xr.DataArray(np.random.rand(180, 360), coords=[lat_old, lon_old], dims=["lat", "lon"])

# Perform conservative regridding to (90, 144)
burned_area_regridded = conservative_regrid_manual(burned_area_data, 90, 144)

# Verify total conservation (this is an approximation)
old_total_ba = burned_area_data.sum()
new_total_ba = burned_area_regridded.sum()

print(f"Total BA before regridding: {old_total_ba.values}")
print(f"Total BA after regridding: {new_total_ba.values}")

# Optional: Save the regridded data to a NetCDF file
burned_area_regridded.to_netcdf("regridded_burned_area.nc")
