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
#include <trackreco/PHRaveVertexing.h>

R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libPHTpcTracker.so)

#endif

#if __cplusplus >= 201703L

#include <G4_Magnet.C>
//#include <trackreco/ActsEvaluator.h>
#include <trackreco/MakeActsGeometry.h>
#include <trackreco/PHActsSourceLinks.h>
#include <trackreco/PHActsTracks.h>
#include <trackreco/PHActsTrkFitter.h>
//#include <trackreco/PHActsTrkProp.h>
////#include <trackreco/PHActsVertexFinder.h>
////#include <trackreco/PHActsVertexFitter.h>
#include <trackreco/PHTpcResiduals.h>

#endif

using namespace std;

int Fun4All_G4_Readback(){

  int verbosity = 4;
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
  map<string, int> reconstructionChannel;
  reconstructionChannel["D02K-pi+"] = 0;
  reconstructionChannel["D02K+pi-"] = 0;
  reconstructionChannel["Bs2Jpsiphi"] = 1;
  reconstructionChannel["Bd2D-pi+"] = 0;
  reconstructionChannel["Upsilon"] = 0;
  reconstructionChannel["testSpace"] = 0;

  const int numberOfActiveRecos = accumulate( begin(reconstructionChannel), end(reconstructionChannel), 0, 
                                              [](const int previous, const pair<const string, int>& element) 
                                              { return previous + element.second; });

  if (numberOfActiveRecos != 1)
  {
    if (numberOfActiveRecos == 0) printf("*\n*\n*\n* You have not enabled any reconstruction, exiting!\n*\n*\n*\n");
    else
    {
        printf("*\n*\n*\n* You have more than one active reconstruction, exiting!\n*\n*\n*\n");
        std::map<string, int>::iterator it;
        for ( it = reconstructionChannel.begin(); it != reconstructionChannel.end(); it++ )
            printf("You have set channel %s to %u\n", it->first.c_str(), it->second);
    }
    delete se;
    gSystem->Exit(0);
    return 1;
  }

  //--------------
  // IO management
  //--------------
  // Hits file
  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  string fileList;
  if (reconstructionChannel["D02K-pi+"] or reconstructionChannel["D02K+pi-"] or reconstructionChannel["testSpace"]) fileList = "fileList_d2kpi.txt";
  if (reconstructionChannel["Bs2Jpsiphi"]) fileList = "fileList_bs2jpsiphi.txt";
  if (reconstructionChannel["Bd2D-pi+"] or reconstructionChannel["Upsilon"]) fileList = "fileList_bbbar.txt";
  //if (reconstructionChannel["testSpace"]) fileList = "fileList_MDC.txt"; 
  hitsin->AddListFile(fileList.c_str());
  se->registerInputManager(hitsin);


  bool use_act_tracking = false;
  if (use_act_tracking)
  {
    #if __cplusplus >= 201703L
    bool SC_CALIBMODE = true;
    /// Geometry must be built before any Acts modules
    MakeActsGeometry* geom = new MakeActsGeometry();
    geom->Verbosity(verbosity);
    //geom->setMagField(G4MAGNET::magfield);
    //geom->setMagFieldRescale(G4MAGNET::magfield_rescale);
    se->registerSubsystem(geom);

    /// Always run PHActsSourceLinks and PHActsTracks first, to convert TrkRClusters and SvtxTracks to the Acts equivalent
    PHActsSourceLinks* sl = new PHActsSourceLinks();
    sl->Verbosity(verbosity);
    sl->setMagField(G4MAGNET::magfield);
    sl->setMagFieldRescale(G4MAGNET::magfield_rescale);
    se->registerSubsystem(sl);

    PHActsTracks* actsTracks = new PHActsTracks();
    actsTracks->Verbosity(verbosity);
    se->registerSubsystem(actsTracks);

    PHActsTrkFitter* actsFit = new PHActsTrkFitter();
    actsFit->Verbosity(verbosity);
    actsFit->doTimeAnalysis(false);
    /// If running with distortions, fit only the silicon+MMs first
    actsFit->fitSiliconMMs(SC_CALIBMODE);
    se->registerSubsystem(actsFit);
  
    if (SC_CALIBMODE)
    {
       /// run tpc residual determination with silicon+MM track fit
       PHTpcResiduals* residuals = new PHTpcResiduals();
       residuals->Verbosity(verbosity);
       se->registerSubsystem(residuals);
    }
   /*
    PHActsVertexFinder* vtxer = new PHActsVertexFinder();
    vtxer->Verbosity(verbosity);
    se->registerSubsystem(vtxer);
    */
    #endif
  }


  bool use_rave_vertexing = false;
  string raveVertexName = "SvtxVertexMapRave";
  if (use_rave_vertexing)
  {
    PHRaveVertexing* rave = new PHRaveVertexing();
    //https://rave.hepforge.org/trac/wiki/RaveMethods
    rave->set_vertexing_method("kalman-smoothing:1");
    rave->set_over_write_svtxvertexmap(false);
    rave->set_svtxvertexmaprefit_node_name(raveVertexName.c_str());
    rave->Verbosity(verbosity);
    se->registerSubsystem(rave);
  }

  //General configurations
  KFParticle_sPHENIX *kfparticle = new KFParticle_sPHENIX();

  const int nEvents = 5e4;

  //Use rave vertexing to construct PV
  if (use_rave_vertexing) kfparticle->setVertexMapNodeName(raveVertexName.c_str());

  kfparticle->setMinimumTrackPT(0.1);
  kfparticle->setMinimumTrackIPchi2(10);
  kfparticle->setMaximumTrackchi2nDOF(1.5);
  kfparticle->setMaximumVertexchi2nDOF(2);
  kfparticle->setMaximumDaughterDCA(0.03);
  kfparticle->setFlightDistancechi2(80);
  kfparticle->setMinDIRA(0.8);
  kfparticle->setMotherPT(0);

  kfparticle->saveDST(1);
  kfparticle->saveOutput(1);
  kfparticle->doTruthMatching(1);
  kfparticle->getDetectorInfo(0);

  std::pair<std::string, int> daughterList[99];
  std::pair<std::string, int> intermediateList[99];
  std::pair<float, float> intermediateMassRange[99];
  int nIntTracks[99];
  float intPt[99];

  if (fileList == "fileList_bbbar.txt") 
    kfparticle->doTruthMatching(0); //I don't think these events have truth variables

  //D2Kpi reco
  if (reconstructionChannel["D02K-pi+"]
  or  reconstructionChannel["D02K+pi-"])
  {
      kfparticle->setMotherName("D0");  
      kfparticle->setMinimumMass(1.7);
      kfparticle->setMaximumMass(2.0);
      kfparticle->setNumberOfTracks(2);
      kfparticle->setMotherIPchi2(50);
    
      kfparticle->constrainToPrimaryVertex(true);
      kfparticle->hasIntermediateStates(false);
      kfparticle->getChargeConjugate(false);

      if (reconstructionChannel["D02K-pi+"])
      {
        daughterList[0] = make_pair("kaon", -1);
        daughterList[1] = make_pair("pion", +1);
        kfparticle->setOutputName("outputData_D02Kmpip_example.root");
      }
      else
      {
        daughterList[0] = make_pair("kaon", +1);
        daughterList[1] = make_pair("pion", -1);
        kfparticle->setOutputName("outputData_D0bar2Kppim_example.root");
      }
  }


  //Bs2Jpsiphi reco
  if (reconstructionChannel["Bs2Jpsiphi"])
  {
      kfparticle->setMotherName("Bs0");  
      kfparticle->setMinimumMass(4.8);
      kfparticle->setMaximumMass(6.0);
      kfparticle->setNumberOfTracks(4);
     
      kfparticle->constrainToPrimaryVertex(true);
      kfparticle->hasIntermediateStates(true);
      kfparticle->constrainIntermediateMasses(true);
      kfparticle->setNumberOfIntermediateStates(2);
    
      intermediateList[0] = make_pair("J/psi", 0);
      daughterList[0]     = make_pair("electron", -1);
      daughterList[1]     = make_pair("electron", +1);
      intermediateMassRange[0] = make_pair(2.5, 3.2);
      nIntTracks[0] = 2;
      intPt[0] = 0;
    
      intermediateList[1] = make_pair("phi", 0);
      daughterList[2]     = make_pair("kaon", -1);
      daughterList[3]     = make_pair("kaon", +1);
      intermediateMassRange[1] = make_pair(0.9, 1.2);
      nIntTracks[1] = 2;
      intPt[1] = 0;

      kfparticle->setOutputName("outputData_Bs2Jpsiphi_example.root");
  }


  //Bd2D-pi+ reco
  if (reconstructionChannel["Bd2D-pi+"])
  {
      kfparticle->setMinimumMass(4.8);
      kfparticle->setMaximumMass(6.0);
      kfparticle->setNumberOfTracks(4);

      kfparticle->hasIntermediateStates(true);
      kfparticle->setNumberOfIntermediateStates(1);

      intermediateList[0] = make_pair("D-", -1);
      daughterList[0]     = make_pair("kaon", +1);
      daughterList[1]     = make_pair("pion", -1);
      daughterList[2]     = make_pair("pion", -1);
      intermediateMassRange[0] = make_pair(1.0, 3.0);
      nIntTracks[0] = 3;
      intPt[0] = 0.;

      daughterList[3] = make_pair("pion", +1);

      kfparticle->setOutputName("outputData_Bd2Dmpip_example.root");
  }


  //Upsilon reco
  if (reconstructionChannel["Upsilon"])
  {
      kfparticle->setMotherName("Upsilon");  
      kfparticle->setMinimumMass(9);
      kfparticle->setMaximumMass(11);
      kfparticle->setNumberOfTracks(8); 
     
      kfparticle->hasIntermediateStates(true);
      kfparticle->setNumberOfIntermediateStates(2);
     
      intermediateList[0] = make_pair("B0", 0);
      daughterList[0]     = make_pair("kaon", +1);
      daughterList[1]     = make_pair("pion", -1);
      daughterList[2]     = make_pair("pion", -1);
      daughterList[3]     = make_pair("pion", +1);
      intermediateMassRange[0] = make_pair(4.8, 6);
      nIntTracks[0] = 4;
      intPt[0] = 0;
     
      intermediateList[1] = make_pair("B0", 0);
      daughterList[4] = daughterList[0];
      daughterList[5] = daughterList[1];
      daughterList[6] = daughterList[2];
      daughterList[7] = daughterList[3];
      intermediateMassRange[1] = make_pair(4.8, 6);
      nIntTracks[1] = 4;
      intPt[1] = 0;

      kfparticle->setOutputName("outputData_Upsilon_example.root");
  }


  //testSpace
  if (reconstructionChannel["testSpace"])
  {
      kfparticle->hasIntermediateStates(false);
      kfparticle->constrainToPrimaryVertex(false);
      kfparticle->getChargeConjugate(false);

      kfparticle->setMotherName("Upsilon");  
      kfparticle->setMinimumMass(0);
      kfparticle->setMaximumMass(11);

      kfparticle->setMinimumTrackIPchi2(0.); //Tracks should point back to the PV
      kfparticle->setMaximumDaughterDCA(0.02);

      kfparticle->setNumberOfTracks(2);
      daughterList[0] = make_pair("electron", +1);
      daughterList[1] = make_pair("electron", -1);

      kfparticle->setOutputName("testSpace.root");
  }


  //More general setup
  kfparticle->setDaughters( daughterList );
  kfparticle->setIntermediateStates( intermediateList );
  kfparticle->setIntermediateMassRange( intermediateMassRange );
  kfparticle->setNumberTracksFromIntermeditateState( nIntTracks );
  kfparticle->setIntermediateMinPT( intPt );

  se->registerSubsystem(kfparticle);
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
