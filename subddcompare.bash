#!/bin/bash -e

# Script to compare daily diagnostic data with monthly data using NetCDF tools (ncks, ncap2).
# Generalized for different diagnostics, species, and number of days per month.
# Script assumes modelE subdaily output with a model time step - 1

# Define diagnostics to loop over (can be modified)
diags="aij"  # e.g., "aijl", "taij", "taijl"

# Define the list of species (variables) to process
species="gustiwind"  # e.g., "z", "aerosol", "temperature"

# If the monthly diag is saved under a different name, define it here, otherwise comment out
monthly_specie="gusti"

# Specify the month and year
year=2010
month=01  # Set the numeric month (e.g., 01 for January, 02 for February, etc.)

# Convert numeric month to a 3-character string for the monthly file
months_str=("JAN" "FEB" "MAR" "APR" "MAY" "JUN" "JUL" "AUG" "SEP" "OCT" "NOV" "DEC")
month_str=${months_str[$(expr $month - 1)]}

# Define modelE output path 
path="/discover/nobackup/projects/giss_ana/users/kmezuman/pyrE/pyrE_standalone/pyrE_inputs/E6F40tE3km/"
simname="E6F40tE3km"

# Dynamically determine the number of days in the specified month and year
days_in_month=$(cal $month $year | awk 'NF {DAYS = $NF}; END {print DAYS}')

# Temporary output file to store combined data
outputfile="out.nc"
echo "Output file path: ${outputfile}"

# Loop over each diagnostic type
for diag in ${diags}; do

  # Remove the output file if it exists from previous runs
  rm -f ${outputfile}

  # Define the monthly file path for comparison (using 3-character month and year)
  monthly="${path}${month_str}${year}.${diag}${simname}.nc"

  # Echo the monthly file being compared (for logging purposes)
  echo "Comparing with monthly file: ${monthly}"

  # Loop over each species
  for specie in ${species}; do

    # Loop over all days in the dynamically determined number of days for the month
    for day in $(seq 1 $days_in_month); do

      # Define the sub-daily file path (using 2-digit month and day format)
      subdd="${path}${year}`printf %02d ${month}``printf %02d ${day}`.${diag}1${simname}.nc"

      # Echo the current file being processed (for logging purposes)
      echo "Processing file: ${subdd}"

      # Use ncks to append the current day's data for the specific species to the output file
      # -h: no header output, -A: append, -v: extract variable
      echo "Appending data from ${subdd} to ${outputfile}"
      ncks -h -A -v ${specie} ${subdd} ${outputfile}

    done

    # Use ncap2 to calculate the mean for the current species across all sub-daily files
    # -s: script to perform the operation, "avg()" calculates the mean, result saved to "subddmean.nc"
    echo "Calculating mean for ${specie} and saving to subddmean.nc"
    ncap2 -h -O -s "mean=avg(${specie})" ${outputfile} "subddmean.nc"

    # Append the monthly data to the file containing the mean of sub-daily data
    echo "Appending monthly data from ${monthly} to subddmean.nc"
    ncks -h -A -v ${monthly_specie} ${monthly} "subddmean.nc"

  done

done
