// $Id: $

/*!
 * \file Fun4All_KFParticle_basic.C
 * \brief 
 * \author Cameron Dean <cdean@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

#include "GlobalVariables.C"
//#include "G4Setup_sPHENIX.C"
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libfun4all.so)

#endif

using namespace std;

int Fun4All_KFParticle_standardContainers(){

  int verbosity = 1;

  //---------------
  // Load libraries
  //---------------
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4dst.so");
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  //---------------
  // Choose reco
  //---------------

  //--------------
  // IO management
  //--------------
  // Hits file
  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  hitsin->AddListFile("fileList_bs2jpsiphi.txt");
  se->registerInputManager(hitsin);

  const int nEvents = 5e4;

  Enable::DSTOUT = false;
  bool save_nTuple = true;
  bool build_jpsi_container = true;
  bool build_d0_container = true;
  bool build_ks0_container = true;
  bool build_phi_container = true;

  DstOut::OutputDir = ".";
  DstOut::OutputFile = "testContainerOutput.root";

  string recoName;
  float minTrackPT = 0.5;
  float minTrackIPchi2 = 15;
  float maxTrackchi2nDOF = 1.5;
  float maxDaughterDCA = 0.02;
  float maxVertexchi2nDOF = 2;
  float minMotherPT = 0;

  /*
   * J/psi Contianer
   */
  //General configurations
  if (build_jpsi_container)
  {
    KFParticle_sPHENIX *kfparticle_jpsi = new KFParticle_sPHENIX();

    kfparticle_jpsi->saveDST(1);
    kfparticle_jpsi->saveOutput(0);
    kfparticle_jpsi->doTruthMatching(0);
    kfparticle_jpsi->getDetectorInfo(0);
    kfparticle_jpsi->hasIntermediateStates(false);
    kfparticle_jpsi->constrainToPrimaryVertex(false);
    kfparticle_jpsi->getChargeConjugate(false);
    kfparticle_jpsi->setVertexMapNodeName("SvtxVertexMap");

    //Track selection
    kfparticle_jpsi->setMinimumTrackPT(minTrackPT);
    kfparticle_jpsi->setMinimumTrackIPchi2(minTrackIPchi2);
    kfparticle_jpsi->setMaximumTrackchi2nDOF(maxTrackchi2nDOF);
    kfparticle_jpsi->setNumberOfTracks(2);
    std::pair<std::string, int> daughterList_jpsi[99];
    daughterList_jpsi[0] = make_pair("electron", +1);
    daughterList_jpsi[1] = make_pair("electron", -1);
    kfparticle_jpsi->setDaughters(daughterList_jpsi);

    //Vertex Selection
    kfparticle_jpsi->setMaximumVertexchi2nDOF(maxVertexchi2nDOF);
    kfparticle_jpsi->setMaximumDaughterDCA(maxDaughterDCA);

    //Mother selection
    recoName = "J/psi";
    kfparticle_jpsi->setMotherName(recoName.c_str());  
    kfparticle_jpsi->setContainerName(recoName.c_str());
    kfparticle_jpsi->setMinimumMass(2.8);
    kfparticle_jpsi->setMaximumMass(3.4);
    kfparticle_jpsi->setMotherPT(minMotherPT);

    if (save_nTuple) kfparticle_jpsi->saveOutput(1);
    if (save_nTuple) kfparticle_jpsi->setOutputName("stdContainer_jpsi.root");
    se->registerSubsystem(kfparticle_jpsi);
  }
  /*
   * D0 Contianer
   */
  if (build_d0_container)
  {
    KFParticle_sPHENIX *kfparticle_d0 = new KFParticle_sPHENIX();

    kfparticle_d0->saveDST(1);
    kfparticle_d0->saveOutput(0);
    kfparticle_d0->doTruthMatching(0);
    kfparticle_d0->getDetectorInfo(0);
    kfparticle_d0->hasIntermediateStates(false);
    kfparticle_d0->constrainToPrimaryVertex(false);
    kfparticle_d0->getChargeConjugate(true);
    kfparticle_d0->setVertexMapNodeName("SvtxVertexMap");

    kfparticle_d0->setMinimumTrackPT(minTrackPT);
    kfparticle_d0->setMinimumTrackIPchi2(minTrackIPchi2);
    kfparticle_d0->setMaximumTrackchi2nDOF(maxTrackchi2nDOF);
    kfparticle_d0->setNumberOfTracks(2);
    std::pair<std::string, int> daughterList_d0[99];
    daughterList_d0[0] = make_pair("kaon", -1);
    daughterList_d0[1] = make_pair("pion", +1);
    kfparticle_d0->setDaughters(daughterList_d0);

    kfparticle_d0->setMaximumVertexchi2nDOF(maxVertexchi2nDOF);
    kfparticle_d0->setMaximumDaughterDCA(maxDaughterDCA);

    recoName = "D0";
    kfparticle_d0->setMotherName(recoName.c_str());  
    kfparticle_d0->setContainerName(recoName.c_str());
    kfparticle_d0->setMinimumMass(1.6);
    kfparticle_d0->setMaximumMass(2.0);

    if (save_nTuple) kfparticle_d0->saveOutput(1);
    if (save_nTuple) kfparticle_d0->setOutputName("stdContainer_d0.root");
    se->registerSubsystem(kfparticle_d0);
  }

  /*
   * Ks0 Contianer
   */
  if (build_ks0_container)
  {
    KFParticle_sPHENIX *kfparticle_ks0 = new KFParticle_sPHENIX();

    kfparticle_ks0->saveDST(1);
    kfparticle_ks0->saveOutput(0);
    kfparticle_ks0->doTruthMatching(0);
    kfparticle_ks0->getDetectorInfo(0);
    kfparticle_ks0->hasIntermediateStates(false);
    kfparticle_ks0->constrainToPrimaryVertex(false);
    kfparticle_ks0->getChargeConjugate(false);
    kfparticle_ks0->setVertexMapNodeName("SvtxVertexMap");

    kfparticle_ks0->setMinimumTrackPT(minTrackPT);
    kfparticle_ks0->setMinimumTrackIPchi2(minTrackIPchi2);
    kfparticle_ks0->setMaximumTrackchi2nDOF(maxTrackchi2nDOF);
    kfparticle_ks0->setNumberOfTracks(2);
    std::pair<std::string, int> daughterList_ks0[99];
    daughterList_ks0[0] = make_pair("pion", -1);
    daughterList_ks0[1] = make_pair("pion", +1);
    kfparticle_ks0->setDaughters(daughterList_ks0);

    kfparticle_ks0->setMaximumVertexchi2nDOF(maxVertexchi2nDOF);
    kfparticle_ks0->setMaximumDaughterDCA(maxDaughterDCA);

    recoName = "KS0";
    kfparticle_ks0->setMotherName(recoName.c_str());
    kfparticle_ks0->setContainerName(recoName.c_str());
    kfparticle_ks0->setMinimumMass(0.4);
    kfparticle_ks0->setMaximumMass(0.6);

    if (save_nTuple) kfparticle_ks0->saveOutput(1);
    if (save_nTuple) kfparticle_ks0->setOutputName("stdContainer_ks0.root");
    se->registerSubsystem(kfparticle_ks0);
  }

  /*
   * phi Contianer
   */
  if (build_phi_container)
  { 
    KFParticle_sPHENIX *kfparticle_phi = new KFParticle_sPHENIX();

    kfparticle_phi->saveDST(1);
    kfparticle_phi->saveOutput(0);
    kfparticle_phi->doTruthMatching(0);
    kfparticle_phi->getDetectorInfo(0);
    kfparticle_phi->hasIntermediateStates(false);
    kfparticle_phi->constrainToPrimaryVertex(false);
    kfparticle_phi->getChargeConjugate(false);
    kfparticle_phi->setVertexMapNodeName("SvtxVertexMap");

    kfparticle_phi->setMinimumTrackPT(minTrackPT);
    kfparticle_phi->setMinimumTrackIPchi2(minTrackIPchi2);
    kfparticle_phi->setMaximumTrackchi2nDOF(maxTrackchi2nDOF);
    kfparticle_phi->setNumberOfTracks(2);
    std::pair<std::string, int> daughterList_phi[99];
    daughterList_phi[0] = make_pair("kaon", -1);
    daughterList_phi[1] = make_pair("kaon", +1);
    kfparticle_phi->setDaughters(daughterList_phi);
 
    kfparticle_phi->setMaximumVertexchi2nDOF(maxVertexchi2nDOF);
    kfparticle_phi->setMaximumDaughterDCA(maxDaughterDCA);

    recoName = "phi";
    kfparticle_phi->setMotherName(recoName.c_str());
    kfparticle_phi->setContainerName(recoName.c_str());
    kfparticle_phi->setMinimumMass(0.9);
    kfparticle_phi->setMaximumMass(1.1);

    if (save_nTuple) kfparticle_phi->saveOutput(1);
    if (save_nTuple) kfparticle_phi->setOutputName("stdContainer_phi.root");
    se->registerSubsystem(kfparticle_phi);
  }

  if (Enable::DSTOUT)
  {
    string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
    //if (Enable::DSTOUT_COMPRESS) DstCompress(out);
    se->registerOutputManager(out);
  }

  //-----------------
  // Event processing
  //-----------------
  if (nEvents < 0)
    return 0;

  se->run(nEvents);

  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
