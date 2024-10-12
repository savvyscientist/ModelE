#!/bin/bash -e
#file names to use:
diags="aijl" # aij taij taijl"
outputfile="out.nc"
for diag in ${diags}; do # loop over all requested diag
  rm ${outputfile} 
  species='z' # a list of species to test
  for specie in ${species}; do # loop ovrer all requested subdd
    for day in {1..31}; do # loop over all days in January
      subdd="/discover/nobackup/projects/giss_ana/users/kmezuman/MASTER/FIRE/AEROCOM/OMAEQSAM_AEROBBs/200601"`printf %02d ${day}`"."${diag}"h48OMAEQSAM_AEROBBs.nc"
      echo ${subdd}
      ncks -h -A -v ${specie} ${subdd} ${outputfile}  #-O:overwrite an existing file -s: String format for text output -v: variable extraction list
    done
    ncap2 -h -O -s "mean=avg(${specie})" ${outputfile} "subddmean.nc"
    monthly="/discover/nobackup/projects/giss_ana/users/kmezuman/MASTER/FIRE/AEROCOM/OMAEQSAM_AEROBBs/JAN2006."${diag}"OMAEQSAM_AEROBBs.nc"
    echo ${monthly}
    ncks -h -A -v ${specie} ${monthly} "subddmean.nc"
  done
done
