// these include guards are not really needed, but if we ever include this
// file somewhere they would be missed and we will have to refurbish all macros
#ifndef MACRO_FUN4ALLG4SLOPECAL_C
#define MACRO_FUN4ALLG4SLOPECAL_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <calib_emc_pi0/CaloCalibEmc_Pi0.h>

R__LOAD_LIBRARY(libfun4all.so)

void Fun4All_G4_Pi0_Tbt(
    const int nEvents = 1,
    const string &inputClustersFileList = "/sphenix/user/jfrantz/caloCalib/xaa",
    const string &outputFile = "test1")
{
  // this convenience library knows all our i/o objects so you don't
  // have to figure out what is in each dst type
  gSystem->Load("libg4dst.so");
  gSystem->Load("libcalibCaloEmc_pi0.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity();  // set it to 1 if you want event printouts

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTClusters");
  in->AddListFile(inputClustersFileList);
  se->registerInputManager(in);

  string outputfile = outputFile + "_g4cemc_eval.root";

  CaloCalibEmc_Pi0 *eval = new CaloCalibEmc_Pi0("CEMC_CALIB_PI0", outputfile);
  //  eval->Verbosity(verbosity);
  se->registerSubsystem(eval);

  se->run(nEvents);
  se->End();
  delete se;
  gSystem->Exit(0);
}

#endif  //MACRO_FUN4ALLG4SLOPECAL_C
