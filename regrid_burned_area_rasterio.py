
import numpy as np
import xarray as xr
import rasterio
from rasterio.warp import reproject, Resampling

def regrid_burned_area_rasterio(burned_area_data, target_lat, target_lon):
    """
    Regrid high-resolution burned area data to a lower-resolution grid using rasterio.

    Parameters:
    - burned_area_data: xarray.DataArray
        The high-resolution burned area data to be regridded.
    - target_lat: 
        The latitude array of the target resolution.
    - target_lon: 
        The longitude array of the target resolution.

    Returns:
    - regridded_data: xarray.DataArray
        The regridded burned area data.
    """

    # Get source dimensions and resolution
    src_transform = burned_area_data.attrs['transform']
    src_crs = burned_area_data.attrs['crs']
    src_height, src_width = burned_area_data.shape
    src_bounds = rasterio.transform.array_bounds(src_height, src_width, src_transform)
    
    # Define the target transform and shape
    target_transform = rasterio.transform.from_bounds(
        west=target_lon.min(),
        south=target_lat.min(),
        east=target_lon.max(),
        north=target_lat.max(),
        width=len(target_lon),
        height=len(target_lat)
    )

    target_shape = (len(target_lat), len(target_lon))

    # Initialize the destination array
    regridded_data_values = np.zeros(target_shape, dtype=burned_area_data.dtype)

    # Perform the reprojection/resampling
    reproject(
        source=burned_area_data.values,
        destination=regridded_data_values,
        src_transform=src_transform,
        src_crs=src_crs,
        src_bounds=src_bounds,
        dst_transform=target_transform,
        dst_crs=src_crs,
        resampling=Resampling.sum
    )

    # Create the xarray.DataArray for the output
    regridded_data = xr.DataArray(
        regridded_data_values,
        coords=[target_lat, target_lon],
        dims=['lat', 'lon'],
        attrs=burned_area_data.attrs
    )

    return regridded_data
