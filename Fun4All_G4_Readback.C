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
    const int nEvents = 5e4,
    const char *inputFile = "/sphenix/data/data03/sphnxpro/user/jinhuang/HF-production-meson-pp200-tracking/D0_piK/D0_piK_3001.cfg.root",
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
  hitsin->AddListFile("fileList_zhaozhong.txt");
  //hitsin->AddListFile("fileList.txt");
  se->registerInputManager(hitsin);

  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX();

  kfparticle->setNumberOfTracks(4);
  std::pair<std::string, int> daughterList[99];
  daughterList[0] = make_pair("electron", -1);
  daughterList[1] = make_pair("electron", +1);
  daughterList[2] = make_pair("kaon", -1);
  daughterList[3] = make_pair("kaon", +1);
  kfparticle->setDaughters( daughterList );

  kfparticle->hasIntermediateStates(true);
  std::pair<std::string, int> intermediateList[99];
  intermediateList[0] = make_pair("J/psi", 0);
  intermediateList[1] = make_pair("phi", 0);
  kfparticle->setIntermediateStates( intermediateList );
  std::pair<float, float> intermediateMassRange[99];
  intermediateMassRange[0] = make_pair(3, 3.2);
  intermediateMassRange[1] = make_pair(0.9, 1.2);
  kfparticle->setIntermediateMassRange( intermediateMassRange );
  int nIntTracks[99];
  nIntTracks[0] = 2;
  nIntTracks[1] = 2;
  kfparticle->setNumberTracksFromIntermeditateState( nIntTracks );
  float intPt[99];
  intPt[0] = 0;
  intPt[1] = 0;
  kfparticle->setIntermediateMinPT( intPt );

  kfparticle->setMinimumMass(5.0);
  kfparticle->setMaximumMass(5.6);
  /*kfparticle->setMinimumTrackPT(0.10);
  kfparticle->setMinimumTrackIPchi2(15);
  kfparticle->setMaximumTrackchi2nDOF(1.5);
  kfparticle->setMaximumVertexchi2nDOF(5);
  kfparticle->setMaximumDaughterDCA(0.3);
  kfparticle->setFlightDistancechi2(60);
  kfparticle->setMinDIRA(0.8);
  kfparticle->setMotherPT(0);*/
  kfparticle->setMinimumTrackPT(0.0);
  kfparticle->setMinimumTrackIPchi2(0);
  kfparticle->setMaximumTrackchi2nDOF(100);
  kfparticle->setMaximumVertexchi2nDOF(100);
  kfparticle->setMaximumDaughterDCA(3);
  kfparticle->setFlightDistancechi2(0);
  kfparticle->setMinDIRA(0.);
  kfparticle->setMotherPT(0);

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
