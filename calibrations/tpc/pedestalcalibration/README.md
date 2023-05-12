Procedure for processing PRDF pedestal run files:

* Make sure you have locally installed the TPCRawDataTree module located at analysis/TPC/DAQ/TPCRawDataTree/
* Edit the `Fun4All_TPC_UnpackPRDF.C` macro to point the outDir variable to whichever directory you'd like to store the output root files in. The macro is located in the `analysis/TPC/DAQ/macros/` directory
* Run the `run_UnpackPRDF.sh` shell script with all files you want to process as an argument, making sure to set the number of events variable in the shell script to be however many events you need. An example of running the script for all TPC sectors in run 10151:

```
source run_UnpackPRDF.sh /sphenix/user/jinhuang/TPC/commissioning/pedestal/TPC_ebdc*_pedestal-00010151-0000.prdf 
``` 

Procedure for creating plots of each sectors dead channels, pedestal mean, and pedestal rms:
* Edit the `TPC_Channel_QA.C` macro to have the `runNumber` string be the name of the file run you want to process and edit the `fileName` to match wherever you have stored the output root files from processing the PRDFs using TPCRawDataTree as above. Also edit the `saveName` variable at the end of the file to point to wherever you'd like the png images to be stored instead of just an `images` directory wherever you're running the macro
* Run `TPC_Channel_QA.C` with `root -l TPC_Channel_QA.C` and you'll generate 24 pngs, one for each sector in the TPC

Procedure for making root file containing channel averaged pedestal info:
* Run the `run_TPCPedestalCalibration.sh` script from the command line with all root files (output from converting the PRDF files) you want to process as an argument. Example for processing all TPC sectors from run 10151:
```
source run_TPCPedestalCalibration.sh /sphenix/user/rosstom/test/testFiles/TPC_ebdc*_pedestal-00010151-0000.prdf_TPCRawDataTree.root
```
* Should produce a new root file for each sector's channels with information about each channel's alive/dead status, pedestal mean, and pedestal rms

If you have any questions, send Thomas Marshall (@rosstom, rosstom@ucla.edu) or Aditya Dash (@adityadash, adityaprasaddash55@gmail.com) a message on mattermost or over email
