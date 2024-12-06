def GFED5_BA(year, config, lons, lats):
    """
    Read and process monthly GFED5 burned area data.
    Calculates natural BA as: Total - Deforestation - Cropland - Peatland
    
    Args:
        year (int): Year to process
        config (dict): Configuration settings
        lons (np.array): Longitude array
        lats (np.array): Latitude array
    
    Returns:
        tuple: (annual_sum_array, total_formatted_string)
    """
    # Initialize arrays
    nlat = len(lats)
    nlon = len(lons)
    ann_sum = np.zeros((nlat, nlon), dtype=float)
    m2toMha = 1E-6  # Convert m² to Mha
    
    # Process each month
    for month in range(1, 13):
        try:
            # Construct filename (BAYYYYMM.nc)
            filename = f"BA{year}{month:02d}.nc"
            filepath = os.path.join(config['dir_obs_ba'], filename)
            
            if not os.path.exists(filepath):
                logging.warning(f"File {filepath} not found. Skipping month.")
                continue
            
            with nc.Dataset(filepath) as f:
                # Read each component
                components = {
                    'Total': None,
                    'Defo': None,
                    'Crop': None,
                    'Peat': None
                }
                
                for comp in components.keys():
                    if comp in f.variables:
                        data = f.variables[comp][:]
                        # Handle missing/fill values if present
                        if hasattr(f.variables[comp], '_FillValue'):
                            fill_value = f.variables[comp]._FillValue
                            data = np.where(data == fill_value, 0, data)
                        # Handle negative values
                        data = np.where(data < 0, 0, data)
                        components[comp] = data
                    else:
                        logging.warning(f"Variable {comp} not found in file for {month}")
                        components[comp] = np.zeros((nlat, nlon))
                
                # Calculate natural BA
                month_ba = (components['Total'] - 
                          components['Defo'] - 
                          components['Crop'] - 
                          components['Peat'])
                
                # Ensure no negative values
                month_ba = np.where(month_ba < 0, 0, month_ba)
                
                # Convert units and add to annual sum
                month_ba *= m2toMha
                ann_sum += month_ba
                
        except Exception as e:
            logging.error(f"Error processing month {month} of {year}: {str(e)}")
            continue
    
    # Calculate total
    total = np.nansum(ann_sum)
    total_formatted = format(total, '.3e')
    
    return ann_sum, total_formatted

# Example usage:
# ba_sum, ba_total = GFED5_BA(2020, config, lons, lats)
