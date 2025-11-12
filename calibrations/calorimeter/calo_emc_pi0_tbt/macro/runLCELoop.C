#include <calib_emc_pi0/CaloCalibEmc_Pi0.h>

#include <TSystem.h>

void runLCELoop(int nevents = -1, const std::string &ifile="treetest_g4cemc_eval.root", const std::string &ofile="treetest1_g4cemc_eval.root")
{
        gSystem->Load("libcalibCaloEmc_pi0.so");
	CaloCalibEmc_Pi0 obj_LCE("name_objLCE", ofile);
	obj_LCE.InitRun(nullptr);
	obj_LCE.Loop(nevents,ifile);
	obj_LCE.End(nullptr);
//	obj_LCE.FittingHistos(); // This method does not exists in CaloCalibEmc_Pi0
}
