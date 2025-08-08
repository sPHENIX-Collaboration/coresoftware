To run this software:
1. build the MVTXMatchingEfficiencyWithShapes
2. folder MVTX_ME_AuAu_template is prepared for the submission:

2a. in run_MDC2reco.sh chaneg the source path - source the previously build library

2b. in run_MDC2reco.sh change the filepath- each list should contain all subsystems DST_HITSETS - i.e. 6 MVTX, 8 INTT, 24 TPC, 1 TPOT

2c. Fun4All_TrackSeeding_25326.C is the main tracking macro, inside the MVTXMatchingEfficiencyWithShapes is called - update accordingly if there any any changes in the Tracking - note, this package requires SearchInINTT()

2d: myMVTXMEAuAu.job is sibmitter that takes two arguments - R=XXXX which represents teh run number which is also the arhument of the folder: MVTX_ME_AuAu_$(R) meaning, template the _template with run number, second argument Q=YYY is the numebr of segments to run over 
