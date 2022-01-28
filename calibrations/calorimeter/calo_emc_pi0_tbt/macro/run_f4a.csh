#!/bin/csh
echo ' '
echo 'START: '`date`
source ~/setup_sphenix.sh
source ~/setldpath.csh sphenix
source /opt/sphenix/core/bin/setup_local.csh $JFINSTALL
echo $LD_LIBRARY_PATH
root.exe -b <<EOF
.x Fun4All_G4_Pi0_Tbt.C($1,"$2","$3")
.q
EOF
echo ' '
echo ' '
echo 'END: '`date`
echo ' '
