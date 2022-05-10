Various macros and files for running tower slope modules in various situations:

Macro/script for running over MDC2 Hits files (Production 4 only) 
-------------------------------------------
run_calo_fromMDC2Hits_towslope_Fun4All_G4_Calo.C
do_run_calo_fromMd2Hits_towslope.csh

Notes: does not need filelist input (e.g. from CreateFileList.pl), see comments inside, takes single 
integer input, assuming filenames of Production 4 MDC2 Hits files.  For other production
runs beginning of macro/filename parsing part my need up dated. Note that per usual MDC2 
and presumably all sphenix default production/sim file access, only the filename is specified
and lustre system/dcache finds the file and decides/provides appropriate "full path" access method.
Valid filenames can be obtained from CreateFileLists.pl
This macro does not need any other files found in this folder to run,
but can be run in with others present to modify behavior.
Runs all three slope method modules one for each calorimeter 
system in central barrel (Emc & Hcal's both)


Alternative G4_Setup / "common" macros for
Applying Decal/Cal through Simple Calo calibration
Databse File API's : 
-------------------------------------------
G4_CEmc_Spacal.C
G4_HcalIn_ref.C  
G4_HcalOut_ref.C

Notes: place any or all of these in running directory to override those in 
$OFFLINE_MAIN/root_macros (originally from sphenix "macros" repository, common/
In order to specify calibrations and decalibrations using RawTowerCalibration/RawTowerDigitizer
as described in Calo Calibrations meetings Jan-March 2022
See comments especially many detailed comments in G4_CEmc_Spacal.C
about how to do this.  Files needed in most cases can be found below:

                       
Example Simple Database Files for Decal/Cal : 
-------------------------------------------
HCALIN_GainsCalib1.12_PosEtaOnly_NegEta1.0.txt 
HCALIN_GainsCalib1.22.txt                        
HCALOUT_GainsCalib1.22.txt
HCALOUT_GainsCalib1.12_PosEtaOnly_NegEta1.0.txt
emc_corr_files/*root *.C

Notes: see inside G4_Setup macros in previous section.
The last emc_corr_files/ subdirectory also contains macros 
(for Emc which uses a root file format as opposed to Hcal 
which uses a text file format) for generating new decal or calibration 
correction patterns for testing

Macros for doing the actual fitting step 
----------------------
do_eta_fit.C
do_eta_fit2.C		                 

The fitting step (timeslice1 --> timeslice 2 --see above)
should generally be after merging of tower energy distribution 
histogram files (the output of the tower slope module over the 
mdc2/production data files).  The last function called in the 
macro LiteCaloEval::Fit_RelativeShift outputs the
correction which can be read as input in place of the "Example
Simple Database" files above.  For CEMC, this is output 
into the root file output containing all the other root
output (various graphs, etc.) into a tree with in that file 
with the correct name ("emc_corr_tree").  For the Hcals,
for which the correction files should be simple text file 
formats, the correction text file is outputted to:
[output_root_filename].HCAL[IN|OUT]_CORR_TXTFILE.txt as 
an additional output file besides the usual output_root_file.


Older/out of date: (from mdc1 running)
----------------------
Fun4All_G4_SlopeCal.C                           
run_f4a.csh


                                do_run_calo_fromMd2Hits_towslope.csh

                                 

