#!/bin/bash -e

# Script to compare daily diagnostic data with monthly data using NetCDF tools (ncks, ncap2).
# Generalized for different diagnostics, species, and number of days per month.

# Define diagnostics to loop over (can be modified)
diags="aijl"  # e.g., "aij", "taij", "taijl"

# Define the list of species (variables) to process
species="z"  # e.g., "z", "aerosol", "temperature"

# Specify the month and year
year=2006
month=01  # Set the numeric month (e.g., 01 for January, 02 for February, etc.)

# Convert numeric month to a 3-character string for the monthly file
months_str=("JAN" "FEB" "MAR" "APR" "MAY" "JUN" "JUL" "AUG" "SEP" "OCT" "NOV" "DEC")
month_str=${months_str[$(expr $month - 1)]}

# Dynamically determine the number of days in the specified month and year
days_in_month=$(cal $month $year | awk 'NF {DAYS = $NF}; END {print DAYS}')

# Output file to store combined data
outputfile="out.nc"

# Loop over each diagnostic type
for diag in ${diags}; do

  # Remove the output file if it exists from previous runs
  rm -f ${outputfile}

  # Define the monthly file path for comparison (using 3-character month and year)
  monthly="/discover/nobackup/projects/giss_ana/users/kmezuman/MASTER/FIRE/AEROCOM/OMAEQSAM_AEROBBs/${month_str}${year}.${diag}OMAEQSAM_AEROBBs.nc"

  # Echo the monthly file being compared (for logging purposes)
  echo "Comparing with monthly file: ${monthly}"

  # Loop over each species
  for specie in ${species}; do

    # Loop over all days in the dynamically determined number of days for the month
    for day in $(seq 1 $days_in_month); do

      # Define the sub-daily file path (using 2-digit month and day format)
      subdd="/discover/nobackup/projects/giss_ana/users/kmezuman/MASTER/FIRE/AEROCOM/OMAEQSAM_AEROBBs/${year}`printf %02d ${month}``printf %02d ${day}`.${diag}h48OMAEQSAM_AEROBBs.nc"

      # Echo the current file being processed (for logging purposes)
      echo "Processing file: ${subdd}"

      # Use ncks to append the current day's data for the specific species to the output file
      # -h: no header output, -A: append, -v: extract variable
      ncks -h -A -v ${specie} ${subdd} ${outputfile}

    done

    # Use ncap2 to calculate the mean for the current species across all sub-daily files
    # -s: script to perform the operation, "avg()" calculates the mean, result saved to "subddmean.nc"
    ncap2 -h -O -s "mean=avg(${specie})" ${outputfile} "subddmean.nc"

    # Append the monthly data to the file containing the mean of sub-daily data
    ncks -h -A -v ${specie} ${monthly} "subddmean.nc"

  done

done
