Collects tower energy distribution for every tower in the Emcal, IHcal, and OHcal, and fits to exponential form in various energy ranges specific to each calo system.

Compiles and runs over MDC1 with 2 oldest macros in the macro/ folder.  
Compiles and runs over MDC2 with obviously named macros/scripts in the macro/ folder.  


Curent status: contains MDC2 updated relative fitting method using data defined histogram shape.  Currently designed to run in "gain tracing" mode where calibration shift in one time slice file is detected from tower energy distributions, quantified, and correction file is generated for, by fitting its shape with same energy distribution for reference time slice. 

See more documentation in macros directory

Mdc2 Initial Commit
-Justin Frantz 5/10/22


OLD LOG: 
For MDC1 It currently runs over DST files (towers only in them)  from MDC1  : e.g. /sphenix/data/data02/sphnxpro/MDC1/sHijing_HepMC/CaloCluster/data/DST_CALO_CLUSTER_sHijing_0_12fm-0000000001-00009.root

TODO:  convert slope fit parameters into calibration factors and commit to files.  Right now it generates an output file with the various tGraphs and Histograms recording the fit results. 

Initial Commit
-Justin Frantz frantz@ohio.edu 5/18/2021
