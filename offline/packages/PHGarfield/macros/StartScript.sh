#! /bin/bash
emin=400
emax=400
ne=1

bmin=1.15
bmax=1.45
nb=50

amin=0.0
amax=0.2
na=50

bnow=$(awk -v bmin="$bmin" \
             -v bmax="$bmax" \
             -v nb="$nb" \
             -v i="$1" \
'BEGIN {
   print bmin + (bmax - bmin) * i / nb
}')

echo $bnow

output="/direct/phenix+u/workarea/hemmick/code.sphenix/tkh/gas/gasfiles/PART_"$1".gas"

echo $output

echo GasModel $emin $emax $ne $bnow $bnow 1 $amin $amax $na $output
GasModel $emin $emax $ne $bnow $bnow 1 $amin $amax $na $output

echo all done
