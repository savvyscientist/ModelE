def process_lightning_density(year, filepath, var_name):
    """
    Process lightning stroke density from monthly data.
    Maps year to correct indices in time dimension (1-144, where 1 is Jan 2013).
    Converts from strokes/km²/day to strokes/km²/year.
    
    Args:
        year (int): Year to process (2013-2020)
        filepath (str): Path to netCDF file
        var_name (str): Name of the density variable in the file
        
    Returns:
        tuple: (annual_sum_array, total_formatted_string)
        Returns (None, None) if year is out of range
    """
    # Check valid year range
    if year < 2013 or year > 2020:
        logging.error(f"Year {year} out of valid range (2013-2020)")
        return None, None
        
    # Calculate time indices for the requested year
    # Index 1 corresponds to Jan 2013
    start_idx = (year - 2013) * 12  # Starting index for requested year
    time_indices = range(start_idx, start_idx + 12)  # 12 months of data
    
    try:
        with nc.Dataset(filepath) as f:
            # Get dimensions
            lats = f.variables['latitude'][:]
            lons = f.variables['longitude'][:]
            
            # Initialize array for annual sum
            annual_density = np.zeros((len(lats), len(lons)), dtype=float)
            
            # Process each month
            for month, time_idx in enumerate(time_indices):
                # Read density data for the month
                density = f.variables[var_name][time_idx,:,:]
                
                # Get number of days in this month
                ndays = calendar.monthrange(year, month+1)[1]
                
                # Convert from per day to total for month
                monthly_total = density * ndays
                
                # Add to annual sum
                annual_density += monthly_total
            
            # Calculate total (optional)
            total = np.nansum(annual_density)
            total_formatted = format(total, '.3e')
            
            return annual_density, total_formatted
            
    except Exception as e:
        logging.error(f"Error processing lightning density data: {str(e)}")
        raise

# Example usage:
# density_sum, total = process_lightning_density(2015, 
#                                              '/path/to/lightning.nc',
#                                              'stroke_density')
