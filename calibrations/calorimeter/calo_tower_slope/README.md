Collects tower energy distribution for every tower in the Emcal, IHcal, and OHcal, and fits to exponential form in various energy ranges specific to each calo system.

Compiles and runs over MDC1 with macro in the macro/ folder.  

It currently runs over DST files (towers only in them)  from MDC1  : e.g. /sphenix/data/data02/sphnxpro/MDC1/sHijing_HepMC/CaloCluster/data/DST_CALO_CLUSTER_sHijing_0_12fm-0000000001-00009.root

TODO:  convert slope fit parameters into calibration factors and commit to files.  Right now it generates an output file with the various tGraphs and Histograms recording the fit results. 

Initial Commit
-Justin Frantz frantz@ohio.edu 5/18/2021
