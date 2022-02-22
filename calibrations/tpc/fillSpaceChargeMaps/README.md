# SetUp
Start here for setting up environment: [Example_of_using_DST_nodes (Wiki)](https://wiki.bnl.gov/sPHENIX/index.php/Example_of_using_DST_nodes#Building%20a%20package)
- Setup local compilation for bash shel:

```
source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=/sphenix/user/shulga/tpc2019_install
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
```
**Do not forget to change the line:** ```export MYINSTALL=/sphenix/user/shulga/tpc2019_install```

<!---
- Creating package files:
```
CreateSubsysRecoModule.pl --all --overwrite fillSpaceChargeMaps
``` 
(*very useful script providing all files needed for new package*, `be carefull with options`)
--->

- Compilation of the package:
```
mkdir build
cd build
../autogen.sh --prefix=$MYINSTALL
make -j 4
make install
```
- reading first file:




# Workflow:
- Files with G4Hits from Hijing events used for the analysis: 
    - G4Hits_sHijing_0_20fm-0000000002-*

- File with  bunchcrossing id and time (ns) assuming 106ns between bunches and 50kHz collision rate: __timestamps_50kHz.txt__. The file is used to mimic the bunchcrossing;

- Running over G4Hits containers in the files is performed with Fun4All environment. Main code is Fun4All_FillChargesMap_300evts*.C, it is run with run_files.sh, which takes as an input the first and last file number:
```
#This will run first 5 files with G4Hits (100 events per file in the MDC2) and create 5 files 
#with histograms:
source macros/run_files_300evts_AA_MDC2.sh 0 5 
```

-  As soon as files are available the histograms are inside;
- To create bunch of bash files and condor job files to start condor jobs scripts are available:
```
#Creating folders to store all the files:

mkdir Out
mkdir Files
mkdir condor_macros

#Creating 1000s job files to run over G4Hits:
scripts/generate_run_files_300evts_AA_MDC2.py
```
**Do not forget to change the path to your repositories:** 

```export MYINSTALL=/sphenix/user/shulga/tpc2019_install```

```/sphenix/user/shulga/Work/...```

- Scripts above will also generate bash files to submit all jobs, *_all bash scripts created above should be provided executable rights before that_*:
```
../run_all_jobs*  
```

# Adding histograms from all files:
The files contain histograms for 100 events each. Full map is ~10000 events. Thus, maps have to be integrated. 
To make files smaller the Sumw2 arrays/matrices for histograms should not be stored.
The tool to provide this functionality:
```
add_histos.py
```
