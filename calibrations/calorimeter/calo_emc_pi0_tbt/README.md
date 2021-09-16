Collects pi0-area invmass by leading tower for every tower in the Emcal.

Compiles and runs over MDC1 with macro in the macro/ folder.  

It currently runs over DST files (towers only in them)  from MDC1  : e.g. /sphenix/data/data02/sphnxpro/MDC1/sHijing_HepMC/CaloCluster/data/DST_CALO_CLUSTER_sHijing_0_12fm-0000000001-00009.root

TODO: At the time of initial commit it is only creating the invmass histos and doing unoptimized fitting on them and saving small cluster tree for applying the derived correction on the next iteration. The code is also functional for subsequent iterations as well using the Loop() function outside of fun4all.  Macros for running over fun4all and later iterations using Loop() are in the macro directory

Initial Commit
-Justin Frantz frantz@ohio.edu 9/9/2021
