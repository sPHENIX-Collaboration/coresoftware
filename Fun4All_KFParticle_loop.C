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
#include <numeric>

R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4dst.so)
#endif

using namespace std;

int verbosity = 4;

std::pair<std::string, int> daughterList[99];

void runReconstruction( string filePath )
{
  const int nEvents = 1e2;
  string inputFileName = filePath.substr(filePath.size() - 20, 20);

  string nTupleName = "D2Kpi_output_nTuples/outputData_D02Kpi_example_";
  nTupleName += inputFileName.substr(inputFileName.size() - 13, 4);
  nTupleName += ".root";

  printf("Processing input file: %s\n", inputFileName.c_str());
  printf("The output nTuple will be %s\n", nTupleName.substr(nTupleName.size()-35, 35).c_str());

  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  hitsin->AddFile(filePath.c_str());
  se->registerInputManager(hitsin);
 
  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX();

  kfparticle->saveOutput(1);
  kfparticle->doTruthMatching(1);
  kfparticle->getDetectorInfo(0);

  //Track and Vertex Selection
  kfparticle->setMinimumTrackPT(0.1);
  kfparticle->setMinimumTrackIPchi2(10);
  kfparticle->setMaximumTrackchi2nDOF(1.5);
  kfparticle->setMaximumVertexchi2nDOF(2);
  kfparticle->setMaximumDaughterDCA(0.03);

  //Mother selection
  kfparticle->setMotherName("D0");  
  kfparticle->setMinimumMass(1.7);
  kfparticle->setMaximumMass(2.0);
  kfparticle->setNumberOfTracks(2);
  kfparticle->setMotherIPchi2(50);
  kfparticle->setMinDIRA(0.8);
  kfparticle->setMotherPT(0);
  kfparticle->setFlightDistancechi2(80);
  
  //Daughter ID and charge  
  daughterList[0] = make_pair("kaon", -1);
  daughterList[1] = make_pair("pion", +1);
  kfparticle->setDaughters( daughterList );

  kfparticle->setOutputName(nTupleName.c_str());
  se->registerSubsystem(kfparticle);

  //-----------------
  // Event processing
  //-----------------
  if (nEvents < 0)
    return 0;
  se->run(nEvents);
  hitsin->fileclose();
  hitsin->ResetFileList();
  se->ResetNodeTree();
  se->unregisterSubsystem(kfparticle);
  cout << "Written output to: " << nTupleName << endl;
  se->EndRun();
  //delete se;
}

std::fstream& GotoLine(fstream& file, int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

int Fun4All_G4_Readback_Loop(int startLine = 1, int stopLine = 3) //max stopLine is 930
{
  //---------------
  // Load libraries
  //---------------
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4dst.so");

  string fileList = "fileList_d2kpi.txt", inputFilePath;

  for (int l = startLine; l < stopLine; ++l)
  {
    fstream infile(fileList);
    GotoLine(infile, l);
    infile>>inputFilePath;
    runReconstruction(inputFilePath);
  }

  gSystem->Exit(0);
  return 0;
}
