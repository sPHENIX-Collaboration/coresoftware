// these include guards are not really needed, but if we ever include this
// file somewhere they would be missed and we will have to refurbish all macros
#ifndef MACRO_RUNHELIXRESIDUALS_C
#define MACRO_RUNHELIXRESIDUALS_C

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <G4_Tracking.C>

#include "helixResiduals.h"

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libmassRecoAnalysis.so)

void runHelixResiduals(
      const int nEvents = 1,
      const string &inputFile = "G4sPHENIX.root",
      const string &outputFile = "residuals.root"
			)
{
  // this convenience library knows all our i/o objects so you don't
  // have to figure out what is in each dst type
  gSystem->Load("libg4dst.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);  // set it to 1 if you want event printouts

  Fun4AllInputManager *inDST = new Fun4AllDstInputManager("inDST");
  inDST -> AddFile(inputFile);
  se -> registerInputManager(inDST);

  Enable::MICROMEGAS = true;
  TrackingInit();

  helixResiduals *eval = new helixResiduals("eval_residuals", outputFile);
  se->registerSubsystem(eval);
  
  se->run(nEvents);
  se->End();
  
  delete se;
  cout << "Analysis Completed" << endl;
  
  gSystem->Exit(0);
}

#endif  //MACRO_RUNHELIXRESIDUALS_C
