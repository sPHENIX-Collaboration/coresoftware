#!/bin/bash

# This shell script looks in the jobs directory and submits every job it finds there.
# It is intended to pair with clone_job_with_new_file.sh, which makes those files.

for file in `ls jobs/`
do
    condor_submit jobs/$file
done
