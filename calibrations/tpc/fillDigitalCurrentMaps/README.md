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
- DST files with TrkrHitSet containing ADCs from Hijing events used for the analysis: 
    - DST_TRKR_HIT_sHijing_0_20fm-0000000006--*

- To get list of the files:
```
CreateFileList.pl --run 6 --type 4 --nopileup DST_TRKR_HIT G4Hits
```

- File with  bunchcrossing id and time (ns) assuming 106ns between bunches and 50kHz collision rate: __timestamps_50kHz.txt__. The file is used to mimic the bunchcrossing;

- Running over containers in the files is performed with Fun4All environment. Main code is Fun4All_FillDCMap.C, it is run with run_files_AA.sh, which takes as an input the first and last file number:
```
#This will run first 5 files with G4Hits (100 events per file in the MDC2) and create files 
#with histograms:
source macros/run_files_AA.sh 1 2 1508071.0
```

-  As soon as files are available the histograms are inside;
- To create bunch of bash files and condor job files to start condor jobs scripts are available:
```
#Creating folders to store all the files:

mkdir Out
mkdir Files
mkdir condor_macros

#Creating 1000s job files to run over G4Hits:
scripts/generate_files_AA.py
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
