// $Id: $

/*!
 * \file Fun4All_G4_Readback.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

//#include <anatutorial/AnaTutorial.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>
#include <phhepmc/Fun4AllHepMCInputManager.h>
#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4dst.so)

#endif

using namespace std;

int Fun4All_G4_Readback(
    const int nEvents = 1e5,
    //const char *inputFile = "G4sPHENIX.root",
    const char *inputFile = "/sphenix/user/mxliu/offline/B-sim/macros/macros/g4simulations/dir-1/G4sPHENIX-5K.root",
    //const char *inputFile = "/sphenix/data/data03/sphnxpro/user/jinhuang/HF-production-meson-pp200-tracking/D0_piK/D0_piK_3001.cfg.root",
    const char *outputFile = "G4sPHENIX_Readback.root")
{
  //---------------
  // Load libraries
  //---------------
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4dst.so");

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(5);

  //--------------
  // IO management
  //--------------
  // Hits file
  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  //hitsin->fileopen(inputFile);
  //hitsin->AddListFile("fileList_bbbar.txt");
  hitsin->AddListFile("fileList.txt");
  se->registerInputManager(hitsin);

  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX();
  kfparticle->setNumberOfTracks(2);
  kfparticle->setMinimumMass(1.8);
  kfparticle->setMaximumMass(1.925);
  kfparticle->setMinimumTrackPT(0.10);
  kfparticle->setMinimumTrackIPchi2(15);
  kfparticle->setMaximumTrackchi2nDOF(1.5);
  kfparticle->setMaximumVertexchi2nDOF(5);
  kfparticle->setMaximumDaughterDCA(0.3);
  kfparticle->setFlightDistancechi2(60);
  kfparticle->setMinDIRA(0.8);
  kfparticle->setMotherPT(0);
  kfparticle->setFirstDaughter("kaon");
  kfparticle->doTruthMatching(0);
  kfparticle->getDetectorInfo(0);
  se->registerSubsystem(kfparticle);

  //-----------------
  // Event processing
  //-----------------
  if (nEvents < 0)
  {
    return 0;
  }

  se->run(nEvents);
  // se->run(0);

  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
