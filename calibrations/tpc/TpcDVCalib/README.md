# TpcDVCalib
sPHENIX TPC drift velocity calibration module

- Compile analysis module
```
mkdir build
cd build
../autogen.sh --prefix=$MYINSTALL
make -j 4
make install
```

- Go to macro directory, submit jobs in condorJob
**Do not forget to change Initialdir, RunNumber, input DST list in Arguments, DVCorrTag which is related with pre-calib drift velocity, and the number of jobs you want to submit**

```
Initialdir     = /sphenix/u/xyu3/hftg01
Executable     = $(Initialdir)/run_data.sh
RunNumber      = 51103
DVCorrTag      = 2
Arguments      = "./runList/trackrunlist/run$(RunNumber)_$INT(Process,%04d).txt ./runList/run$(RunNumber)_calo.list $(DVCorrTag)"
Queue 9
```

- Output root file saved in Reconstructed/[runnumber]

- Do offline analysis to calibrate drift velocity
```
./offline_dv_calib.sh [runnumber]
```

- Figure saved in figure/ directory
