def modelE_diag(diag, year, config, lons, lats, axyp):
    """
    Calculate annual sum of ModelE diagnostics from monthly files.
    For BA, sums BA_tree/shrub/grass components for each month.
    For non-BA diagnostics, applies area weighting using axyp.
    
    Args:
        diag (str): Diagnostic to process ('BA', 'fireCount', 'CtoG', 'flammability')
        year (int): Year to process
        config (dict): Configuration settings
        lons (np.array): Longitude array
        lats (np.array): Latitude array
        axyp (np.array): Grid cell areas for area weighting
    
    Returns:
        tuple: (annual_sum_array, total_formatted_string)
    """
    months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 
              'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    m2toMha = 1E-6  # Convert m² to Mha
    
    # Initialize arrays
    nlat = len(lats)
    nlon = len(lons)
    ann_sum = np.zeros((nlat, nlon), dtype=float)

    for month in months:
        try:
            monthly_filename = f"{month}{year}.aijnk_CCycle_E6obioF40.nc"
            filepath = os.path.join(config['dir_sim'], monthly_filename)
            
            if not os.path.exists(filepath):
                logging.warning(f"File {filepath} not found. Skipping month.")
                continue
                
            with nc.Dataset(filepath) as f:
                if diag == 'BA':
                    # Sum the three BA components for each month
                    ba_types = ['BA_tree', 'BA_shrub', 'BA_grass']
                    month_sum = np.zeros((nlat, nlon), dtype=float)
                    for ba_type in ba_types:
                        if ba_type in f.variables:
                            ba_data = f.variables[ba_type][:]
                            ba_data *= m2toMha  # Convert to Mha
                            month_sum += ba_data
                        else:
                            logging.warning(f"Variable {ba_type} not found in file for {month}")
                    ann_sum += month_sum
                else:
                    # Handle other diagnostics
                    if diag in f.variables:
                        diag_data = f.variables[diag][:]
                        units = f.variables[diag].units
                        scaling_factor, unit = extract_scaling_factor(units)
                        diag_data *= float(scaling_factor)
                        # Add monthly value to annual sum
                        ann_sum += diag_data
                    else:
                        logging.warning(f"Variable {diag} not found in file for {month}")
        
        except Exception as e:
            logging.error(f"Error processing {month} {year}: {str(e)}")
            continue
    
    # Calculate total - different handling for BA vs other diagnostics
    if diag == 'BA':
        total = np.nansum(ann_sum)  # Simple sum for BA
    else:
        # Apply area weighting for non-BA diagnostics
        total = np.nansum(ann_sum * axyp)  # Area-weighted sum
        
    total_formatted = format(total, '.3e')
    
    return ann_sum, total_formatted

# Example usage:
# For BA:
# ba_sum, ba_total = modelE_diag('BA', year, config, lons, lats, axyp)
# 
# For other diagnostics (with area weighting):
# fire_count_sum, fire_count_total = modelE_diag('fireCount', year, config, lons, lats, axyp)
