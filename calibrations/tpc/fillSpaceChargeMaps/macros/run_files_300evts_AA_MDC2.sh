#!/usr/bin/bash
start=${1?Error: no start file \# given}
stop=${2?Error: no stop file \# given}
bX=${3?Error: no beamX \# given}
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export MYINSTALL=/sphenix/user/shulga/tpc2019_install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
 
for (( f=$start; f<$stop; f++ ));
do 
    let Xstart=f XevtStart=f*100;
    A=$( printf '%05d' $Xstart )
    #B=$( printf '%06d' $Xend )
    #fname="G4Hits_sHijing_0_20fm-0000000002-"$A".root" ;
    fname="G4Hits_sHijing_0_20fm-0000000040-"$A".root" ;
    foutputname="./Files/mdc2_ADCBins_UseFieldMaps_hist_G4Hits_sHijing_0-12fm_bX"$bX"_"$A".root" ;
    #foutputname="./Files/mdc2_ADCBins_NoFieldMaps_hist_G4Hits_sHijing_0-12fm_bX"$bX"_"$A".root" ;
    echo $fname ;
    echo $foutputname ;
    root -l -b -q ./macros/Fun4All_FillChargesMap_300evts_MDC2.C\(1,$XevtStart,$bX,\"$fname\",\"$foutputname\"\)
done

echo all done
