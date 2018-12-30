#! /bin/tcsh -f 

echo "Usage: $0 new_QA_file [ reference_QA_file ]";

set q = '"';

if ($# < 1) then
 echo "Missing parameters: $*"
 exit 1
endif

set new_QA_file = "$q$1$q";
set reference_QA_file = 'NULL';


if ($# >= 2) then
	set reference_QA_file = "$q$2$q";
endif

echo "$0 - New QA file: $new_QA_file";
echo "$0 - Reference QA file: $reference_QA_file";

set macros = (\
	QA_Draw_CEMC_G4Hit.C \
	QA_Draw_CEMC_TowerCluster.C \
	QA_Draw_HCALIN_G4Hit.C \
	QA_Draw_HCALIN_TowerCluster.C \
	QA_Draw_HCALOUT_G4Hit.C \
	QA_Draw_HCALOUT_TowerCluster.C \
	QA_Draw_Sum_Cluster.C \
	QA_Draw_Calorimeter_Sum_TrackProj.C \
	QA_Draw_Calorimeter_Sum_TrackProjEP.C \
);

# imake nstall-data
foreach macro ($macros)
	root -b -q "$macro($new_QA_file, $reference_QA_file)"
end

# last all jet stuff
root -b -q "QA_Draw_Jet_Spectrum.C(${q}h_QAG4SimJet_AntiKt_Truth_r07${q}, $new_QA_file, $reference_QA_file)"
root -b -q "QA_Draw_Jet_Spectrum.C(${q}h_QAG4SimJet_AntiKt_Tower_r07${q}, $new_QA_file, $reference_QA_file)"
root -b -q "QA_Draw_Jet_Summary.C(${q}AntiKt_Tower${q},$new_QA_file, $reference_QA_file)"

echo "$0 - Output plots:";
ls -lh $1*.png;

