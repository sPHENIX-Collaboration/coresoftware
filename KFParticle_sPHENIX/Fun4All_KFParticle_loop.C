// $Id: $

/*!
 * \file Fun4All_KFParticle_loop.C
 * \brief 
 * \author Cameron Dean <cdean@bnl.gov>
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

int verbosity = 0;

std::pair<std::string, int> daughterList[99];

Fun4AllServer* server()
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(verbosity);
  return se;
}

Fun4AllInputManager* inputManager(string filePath)
{
  string inputFileName = filePath.substr(filePath.size() - 20, 20);
  string fileNumber = inputFileName.substr(inputFileName.size() - 13, 4);
  string inputManagerName = "DSTin_" + fileNumber;

  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager(inputManagerName.c_str());
  hitsin->AddFile(filePath.c_str());

  return hitsin; 
}

KFParticle_sPHENIX* eventReconstruction(string filePath)
{
  string inputFileName = filePath.substr(filePath.size() - 20, 20);

  string fileNumber = inputFileName.substr(inputFileName.size() - 13, 4);
  string nTupleName = "D2Kpi_output_nTuples/outputData_D02Kpi_example_";
  nTupleName += fileNumber;
  nTupleName += ".root";
  string kfparticleName = "KFPARTICLE_" + fileNumber;
  printf("The output nTuple will be %s\n", nTupleName.substr(nTupleName.size()-35, 35).c_str());

  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX(kfparticleName.c_str());

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

  return kfparticle;
}

std::fstream& GotoLine(fstream& file, int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}

int Fun4All_KFParticle_loop(int startLine = 2, int stopLine = 5) //max stopLine is 930
{
  //---------------
  // Load libraries
  //---------------
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4dst.so");

  string fileList = "fileList_d2kpi.txt", inputFilePath;

  unsigned int nFiles = stopLine - startLine;
  Fun4AllServer* se[nFiles];
  Fun4AllInputManager* input[nFiles];
  KFParticle_sPHENIX* d2kpi_reco[nFiles];

  const int nEvents = 1e2;

  for (int i = startLine; i < stopLine; ++i)
  {
    fstream infile(fileList);
    GotoLine(infile, i);
    infile>>inputFilePath;

    se[i] = server();

    input[i] = inputManager(inputFilePath);
    se[i]->registerInputManager(input[i]);

    d2kpi_reco[i] = eventReconstruction(inputFilePath);
    se[i]->registerSubsystem(d2kpi_reco[i]);

    //-----------------
    // Event processing
    //-----------------
    if (nEvents < 0)
      return 0;
    se[i]->run(nEvents);

    se[i]->Reset();
  }

  gSystem->Exit(0);
  return 0;
}
