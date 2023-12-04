#!/bin/bash

# This shell script follows the pattern of the valid sample job, with some hints
# about what parts are what, and generates a job by bXing number for all other
# source files in the provided filelist.
# If you are changing the format of the filename, you may need to alter the code
# so that it continues to extract and replace the correct items.
#
# in: (hardcoded) example file, hints about what to replace, and filelist
# out: set of job_NNNNN.job files, each of which can be submitted to condor.

#specify what parts of the sample to replace:
sample_job=condor_run_testfile.job
sample_job_sourcefile=/sphenix/user/shulga/Work/IBF/DistortionMap/Files/Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX10556072.root
sample_job_outfile=Summary_hist_mdc2_UseFieldMaps_AA_event_0_bX10556072
sample_job_outdir=/sphenix/user/rcorliss/distortion_maps/2022.01
sample_job_crossing=10556072

#specify what to replace them with (FFFFF means 'full source path'.  NNNNN means 'bxing number')
this_job_sourcefile=FFFFF
this_job_outfile=Summary_hist_mdc2_UseFieldMaps_AA_event_0_bXNNNNN
this_job_outdir=/sphenix/u/czhang4/distortion_maps/2023.02

#get a list of the files:
filelist=`ls /sphenix/user/shulga/Work/workfest2021_pull/coresoftware/calibrations/tpc/fillSpaceChargeMaps/Files/Summary_hist_mdc2_UseFieldMaps_AA_event_*`


#loop over the elements in that list

for sourcefile in $filelist
do
    #find the root number of the crossing
    crossing=`echo $sourcefile | grep -Eo '[0-9]{6,9}' `
    #make the name of the specific job
    jobname=jobs/job_${crossing}.job
    cp $sample_job $jobname
    echo assigning job: $jobname with crossing: $crossing 
    echo $sample_job_sourcefile $this_job_sourcefile
    #filter the instances of the original filename with the new filename in the specific job
    sed -i "s|${sample_job_sourcefile}|${this_job_sourcefile}|g" $jobname
    sed -i "s|${sample_job_outfile}|${this_job_outfile}|g" $jobname
    sed -i "s|${sample_job_outdir}|${this_job_outdir}|g" $jobname
    sed -i "s|${sample_job_crossing}|${crossing}|g" $jobname
    sed -i "s|FFFFF|${sourcefile}|g" $jobname
    sed -i "s|NNNNN|${crossing}|g" $jobname

    #uncomment to run just one, for testing:
    #exit
done
