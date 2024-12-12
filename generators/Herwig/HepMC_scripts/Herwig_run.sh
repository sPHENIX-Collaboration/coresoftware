#! /bin/bash
filenumber={$3-0}
filetag="_%06d" $filenumber
source /cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/bin/sphenix_setup.sh -n ana
 
/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/Herwig/bin/Herwig run $1 -N $2 -d1 -t $filetag
