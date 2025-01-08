#!/usr/bin/bash
#IN_FILE="/sphenix/user/shulga/Work/TpcPadPlane_v3_coresoftware/macros/detectors/sPHENIX/Files/G4sPHENIX_NewTestGeo_NewRLimits_1evt.root"
#OUT_FILE="/sphenix/user/shulga/Work/TpcPadPlane_v3_coresoftware/macros/detectors/sPHENIX/Files/DC_G4sPHENIX_NewTestGeo_NewRLimits_1evt.root"
#IN_FILE="/sphenix/user/shulga/Work/TpcPadPlane_v3_coresoftware/macros/detectors/sPHENIX/Files/G4sPHENIX_NewTestGeo_1evt.root"
#OUT_FILE="/sphenix/user/shulga/Work/TpcPadPlane_v3_coresoftware/macros/detectors/sPHENIX/Files/DC_G4sPHENIX_NewTestGeo_1evt.root"
#IN_FILE="/sphenix/user/shulga/Work/TpcPadPlane_v3_coresoftware/macros/detectors/sPHENIX/Files/G4sPHENIX_Default_1evts.root"
#OUT_FILE="/sphenix/user/shulga/Work/TpcPadPlane_v3_coresoftware/macros/detectors/sPHENIX/Files/DC_G4sPHENIX_Default_1evt.root"
IN_FILE="DST_TRKR_HIT_sHijing_0_20fm-0000000006-000001.root"
OUT_FILE="/sphenix/user/shulga/Work/IBF/readDigitalCurrents/Files/DC_G4sPHENIX_Default_1evt.root"
root -l -b -q ./macros/Fun4All_FillDCMap.C\(1000,0,1508071,\"$IN_FILE\",\"$OUT_FILE\"\)