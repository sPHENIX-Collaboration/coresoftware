#!/usr/bin/bash
#inputname=${1?Error: no input file given}
#outputbase=${2?Error: no output basename given}
gainname="not_using_gain"

inputname=/sphenix/user/shulga/Work/IBF/DistortionMap/Files/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX10556072.root
outputbase=/direct/star+u/rcorliss/sphenix/generator_output/zerospacecharge_shift_z

source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=/star/u/rcorliss/install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
 
for (( z=-2; z<2; z++ ));
do 
    outputname=${outputbase}_${z} ;
    echo Processing $inputname index $n to output: $outputname;
    echo $foutputname ;
    #signature of the function call is: void generate_distortion_map(const char *inputname, const char* gainName, const char *outputname, const char *ibfName, const char *primName, bool hasSpacecharge=true, bool isAdc=false, int nSteps=500, bool scanSteps=false, float zshift=0){

    root -b -q ./generate_distortion_map.C\(\"${inputname}\",\"$gainname\",\"$outputname\",\"_h_SC_ibf_0\",\"_h_SC_prim_0\",0,0,300,0,${z}\)
done

echo all done
