#!/usr/bin/bash
if [ "$1" == "" ] || [ "$2" == "" ]; then
    echo "Usage: bash_for_condor_submission.sh <inputfile> <outputbasename>"
    echo "or bash_for_condor_submission.sh no"
    exit 1
fi
if [ "$1" == "no" ]; then
echo "Running in quick test mode"
inputname=/sphenix/user/shulga/Work/IBF/DistortionMap/Files/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX10556072.root
outputbase=/direct/star+u/rcorliss/sphenix/generator_output/quick_test
gainname="not_using_gain"
else
inputname=${1?Error: no input file given}
outputbase=${2?Error: no output basename given}
gainname="not_using_gain"
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=/sphenix/u/czhang4/install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
 
for (( n=0; n<30; n++ ));
do 
    outputname=${outputbase}_${n} ;
    echo Processing $inputname index $n to output: $outputname;
    echo $foutputname ;
    root -b -q ./generate_distortion_map.C\(\"${inputname}\",\"$gainname\",\"$outputname\",\"_h_SC_ibf_${n}\",\"_h_SC_prim_${n}\",1\)
done

echo all done
