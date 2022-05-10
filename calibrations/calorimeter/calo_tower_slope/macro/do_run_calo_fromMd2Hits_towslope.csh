#!/bin/csh
echo ' '
echo 'START: '`date`
#setenv HOME /phenix/u/$LOGNAME
#source /etc/csh.login
#source $HOME/.login
#source /opt/phenix/bin/phenix_setup.c 


#@ z = $1 + 500
@ z = $1 + 0

echo hi justin doing job num $z

source /opt/sphenix/core/bin/sphenix_setup.csh -n new

setenv JFINSTALL /phenix/scratch/jfrantz/install/sphenix_bt_mdc2
ls $JFINSTALL/include

source /opt/sphenix/core/bin/setup_local.csh $JFINSTALL

#//setenv JFINSTALL /sphenix/user/jfrantz/newgeom/install
#//setenv LD_LIBRARY_PATH $MYINSTALL/lib:$LD_LIBRARY_PATH
#//set PATH = ($MYINSTALL/bin $path)
#source ~/setldpath.csh . 
#//echo $LD_LIBRARY_PATH
#//echo $MYINSTALL 
echo $PATH

#setenv CALIBRATIONROOT $JFINSTALL/share/calibrations/
#setenv CALIBRATIONROOT /phenix/scratch/jfrantz/install/sphenix_ynewy/share/calibrations/

#cd /direct/phenix+u/apun/sphenixabi/macros/macros/g4simulations_embed
#cd /gpfs/mnt/gpfs04/sphenix/user/apun/gitfork/macros/macros/g4simulations
#cd /direct/phenix+u/apun/sphenixabi/macros/macros/g4simulations
#cd /direct/phenix+u/apun/Island/macros/g4simulations
cd /sphenix/user/jfrantz/caloCalib/framew/calotest/

 root.exe -b <<EOF
   .L run_calo_fromMDC2Hits_towslope_Fun4All_G4_Calo.C
   run_calo_fromMDC2Hits_towslope_Fun4All_G4_Calo(80, $z,"/sphenix/user/jfrantz/caloCalib/framew/condor/v12/out$z.root");
		
.q
EOF
echo ' '
echo ' '
echo 'END: '`date`
echo ' '

#//	"DST_CALO_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000004-$1.root",
#//	"DST_VERTEX_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000004-$1.root",

#  Fun4All_G4_sPHENIX(50,"/sphenix/sim/sim01/production/2016-07-21/single_particle/spacal2d/fieldmap/G4Hits_sPHENIX_pi0_eta0_8GeV-0000.root","/gpfs/mnt/gpfs04/sphenix/user/apun/condor/reduced/reduced_RCBv1_gamma_eta0_16GeV_1evt_$1","","/gpfs/mnt/gpfs04/sphenix/user/apun/condor/withv1/test7/$1.root")

#.L Fun4All_G4_sPHENIX.C
#Fun4All_G4_sPHENIX(50,"/sphenix/sim/sim01/production/2016-07-21/single_particle/spacal2d/fieldmap/G4Hits_sPHENIX_pi0_eta0_8GeV-0000.root","reduced_RCBv1_gamma_eta0_16GeV_1evt_$1","/gpfs/mnt/gpfs02/sphenix/user/apun/condor/withv1/test1/$1.root")

#.L Fun4All_G4_sPHENIX_embed.C
#Fun4All_G4_sPHENIX_embed(30,"/sphenix/sim/sim01/production/2016-07-21/single_particle/spacal2d/fieldmap/G4Hits_sPHENIX_pi0_eta0_8GeV-0000.root","/gpfs/mnt/gpfs02/sphenix/user/apun/condor/reduced/reduced_RCBv1_gamma_eta0_16GeV_1evt_$1","/gpfs/mnt/gpfs02/sphenix/sim/sim01/production/sHijing/2016-12-21/fm_0-4/G4Hits_sHijing_0-4fm_02500_02600.root","/gpfs/mnt/gpfs02/sphenix/user/apun/condor/embedd/withv1/test1/$1.root")





#.L Fun4All_G4_sPHENIX_IA.C
#Fun4All_G4_sPHENIX_IA(50,8.,0.1,"Gamma",1,"/sphenix/sim/sim01/production/2016-07-21/single_particle/spacal2d/fieldmap/G4Hits_sPHENIX_pi0_eta0_8GeV-0000.root","/gpfs/mnt/gpfs04/sphenix/user/apun/condor/reduced/reduced_RCBv1_gamma_eta0_16GeV_1evt_$1","/gpfs/mnt/gpfs04/sphenix/user/apun/condor/island/test4/$1.root")
