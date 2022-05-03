#ifndef MACRO_FUN4ALLG4CALO_C
#define MACRO_FUN4ALLG4CALO_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup_sPHENIX.C>
#include <G4_Bbc.C>
#include <G4_CaloTrigger.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>
#include <G4_Tracking.C>
#include <G4_User.C>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <litecaloeval/LiteCaloEval.h>


#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>


R__LOAD_LIBRARY(libfun4all.so)

// For HepMC Hijing
// try inputFile = /sphenix/sim/sim01/sphnxpro/sHijing_HepMC/sHijing_0-12fm.dat

int Fun4All_wtowSlope_G4_Calo(
    const int nEvents = 10,
      const string &inputFile0 = "/gpfs/mnt/gpfs02/sphenix/sim/sim01/sphnxpro/MDC2/sHijing_HepMC/fm_0_20/pileup/DST_CALO_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000003-00001.root",  //"DST_CALO_G4HIT_sHijing_0_12fm-0000000001-00000.root",
      const string &inputFile1 = "/gpfs/mnt/gpfs02/sphenix/sim/sim01/sphnxpro/MDC2/sHijing_HepMC/fm_0_20/pileup/DST_VERTEX_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000003-00001.root", //"DST_VERTEX_sHijing_0_12fm-0000000001-00000.root",
      const string &outputFile = "G4sPHENIX_calo2",
      const string &outdir = ".")
{
   Fun4AllServer *se = Fun4AllServer::instance();
   se->Verbosity(0);

  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  PHRandomSeed::Verbosity(1);

  /*
  long mtime = gSystem->Now();   
  TString fileOut = outputFile.c_str();
  fileOut += (mtime / 100000 - 8542922)/3;
  string outputFile2 = fileOut.Data();
  */

  string outputFile2 = outputFile.c_str();
  outputFile2 = outputFile2 + ".root";


  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as seed
  // You can either set this to a random value using PHRandomSeed()
  // which will make all seeds identical (not sure what the point of
  // this would be:
  //  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
  // or set it to a fixed value so you can debug your code
  //  rc->set_IntFlag("RANDOMSEED", 12345);

  //===============
  // Input options
  //===============
  // verbosity setting (applies to all input managers)
  Input::VERBOSITY = 1;
  // First enable the input generators
  // Either:
  // read previously generated g4-hits files, in this case it opens a DST and skips
  // the simulations step completely. The G4Setup macro is only loaded to get information
  // about the number of layers used for the cell reco code
  Input::READHITS = true;
  INPUTREADHITS::filename[0] = inputFile0;
  INPUTREADHITS::filename[1] = inputFile1;

  //-----------------
  // Initialize the selected Input/Event generation
  //-----------------
  // This creates the input generator(s)
  InputInit();

  // register all input generators with Fun4All
  InputRegister();

  // set up production relatedstuff
   Enable::PRODUCTION = true;

  //======================
  // Write the DST
  //======================

  Enable::DSTOUT = true;
  Enable::DSTOUT_COMPRESS = false;
  DstOut::OutputDir = outdir;
//  DstOut::OutputFile = outputFile;
  DstOut::OutputFile = outputFile2;


  //Option to convert DST to human command readable TTree for quick poke around the outputs
  //  Enable::DSTREADER = true;

  // turn the display on (default off)
  Enable::DISPLAY = false;

  //======================
  // What to run
  //======================
  // Global options (enabled for all enables subsystems - if implemented)
  //  Enable::ABSORBER = true;
  //  Enable::OVERLAPCHECK = true;
  //  Enable::VERBOSITY = 1;

  Enable::CEMC = true;
  Enable::CEMC_CELL = Enable::CEMC && true;
  Enable::CEMC_TOWER = Enable::CEMC_CELL && true;
  Enable::CEMC_CLUSTER = Enable::CEMC_TOWER && true;
//  Enable::CEMC_EVAL = Enable::CEMC_CLUSTER && true;

  Enable::HCALIN = true;
  Enable::HCALIN_CELL = Enable::HCALIN && true;
  Enable::HCALIN_TOWER = Enable::HCALIN_CELL && true;
  //  Enable::HCALIN_CLUSTER = Enable::HCALIN_TOWER && true;
//  Enable::HCALIN_EVAL = Enable::HCALIN_CLUSTER && true;

  Enable::HCALOUT = true;
  Enable::HCALOUT_CELL = Enable::HCALOUT && true;
  Enable::HCALOUT_TOWER = Enable::HCALOUT_CELL && true;
  // Enable::HCALOUT_CLUSTER = Enable::HCALOUT_TOWER && true;
//  Enable::HCALOUT_EVAL = Enable::HCALOUT_CLUSTER && true;


//  Enable::EPD = false;

//  Enable::GLOBAL_RECO = true;
  //  Enable::GLOBAL_FASTSIM = true;

//  Enable::CALOTRIGGER = Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER && false;

//  Enable::JETS = true;
//  Enable::JETS_EVAL = Enable::JETS && true;

  // HI Jet Reco for p+Au / Au+Au collisions (default is false for
  // single particle / p+p-only simulations, or for p+Au / Au+Au
  // simulations which don't particularly care about jets)
//  Enable::HIJETS = false && Enable::JETS && Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER;

  // 3-D topoCluster reconstruction, potentially in all calorimeter layers
//  Enable::TOPOCLUSTER = true && Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER;
  // particle flow jet reconstruction - needs topoClusters!
//  Enable::PARTICLEFLOW = true && Enable::TOPOCLUSTER;

  //---------------
  // Magnet Settings
  //---------------

  //  const string magfield = "1.5"; // alternatively to specify a constant magnetic field, give a float number, which will be translated to solenoidal field in T, if string use as fieldmap name (including path)
  //  G4MAGNET::magfield = string(getenv("CALIBRATIONROOT")) + string("/Field/Map/sPHENIX.2d.root");  // default map from the calibration database
//  G4MAGNET::magfield_rescale = 1;  // make consistent with expected Babar field strength of 1.4T

  //---------------
  // Pythia Decayer
  //---------------
  // list of decay types in
  // $OFFLINE_MAIN/include/g4decayer/EDecayType.hh
  // default is All:
  // G4P6DECAYER::decayType = EDecayType::kAll;

  // Initialize the selected subsystems
  G4Init();

  //---------------------
  // GEANT4 Detector description
  //---------------------
  if (!Input::READHITS)
  {
    G4Setup();
  }

  //------------------
  // Detector Division
  //------------------

  if (Enable::BBC || Enable::BBCFAKE) Bbc_Reco();

  if (Enable::MVTX_CELL) Mvtx_Cells();
  if (Enable::INTT_CELL) Intt_Cells();
  if (Enable::TPC_CELL) TPC_Cells();
  if (Enable::MICROMEGAS_CELL) Micromegas_Cells();

  if (Enable::CEMC_CELL) CEMC_Cells();

  if (Enable::HCALIN_CELL) HCALInner_Cells();

  if (Enable::HCALOUT_CELL) HCALOuter_Cells();

  //-----------------------------
  // CEMC towering and clustering
  //-----------------------------

  if (Enable::CEMC_TOWER) CEMC_Towers();
  if (Enable::CEMC_CLUSTER) CEMC_Clusters();

  //-----------------------------
  // HCAL towering and clustering
  //-----------------------------

  if (Enable::HCALIN_TOWER) HCALInner_Towers();
  if (Enable::HCALIN_CLUSTER) HCALInner_Clusters();

  if (Enable::HCALOUT_TOWER) HCALOuter_Towers();
  if (Enable::HCALOUT_CLUSTER) HCALOuter_Clusters();

  // if enabled, do topoClustering early, upstream of any possible jet reconstruction
  if (Enable::TOPOCLUSTER) TopoClusterReco();

  if (Enable::DSTOUT_COMPRESS) ShowerCompress();

  //--------------
  // SVTX tracking
  //--------------
  if (Enable::MVTX_CLUSTER) Mvtx_Clustering();
  if (Enable::INTT_CLUSTER) Intt_Clustering();
  if (Enable::TPC_CLUSTER) TPC_Clustering();
  if (Enable::MICROMEGAS_CLUSTER) Micromegas_Clustering();

  if (Enable::TRACKING_TRACK)
  {
    TrackingInit();
    Tracking_Reco();
  }
  //-----------------
  // Global Vertexing
  //-----------------

  if (Enable::GLOBAL_RECO && Enable::GLOBAL_FASTSIM)
  {
    cout << "You can only enable Enable::GLOBAL_RECO or Enable::GLOBAL_FASTSIM, not both" << endl;
    gSystem->Exit(1);
  }
  if (Enable::GLOBAL_RECO)
  {
    Global_Reco();
  }
  else if (Enable::GLOBAL_FASTSIM)
  {
    Global_FastSim();
  }

  //-----------------
  // Calo Trigger Simulation
  //-----------------

  if (Enable::CALOTRIGGER)
  {
    CaloTrigger_Sim();
  }

  //---------
  // Jet reco
  //---------

  if (Enable::JETS) Jet_Reco();
  if (Enable::HIJETS) HIJetReco();

  if (Enable::PARTICLEFLOW) ParticleFlow();

  //----------------------
  // Simulation evaluation
  //----------------------
  string outputroot = outputFile;
  string remove_this = ".root";
  size_t pos = outputroot.find(remove_this);
  if (pos != string::npos)
  {
    outputroot.erase(pos, remove_this.length());
  }

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  if (Enable::PRODUCTION)
  {
    Production_CreateOutputDir();
  }

  if (Enable::DSTOUT)
  {
    string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
    out->AddNode("Sync");
    out->AddNode("EventHeader");
    out->AddNode("TOWER_SIM_HCALIN");
    out->AddNode("TOWER_RAW_HCALIN");
    out->AddNode("TOWER_CALIB_HCALIN");
    out->AddNode("CLUSTER_HCALIN");
    out->AddNode("TOWER_SIM_HCALOUT");
    out->AddNode("TOWER_RAW_HCALOUT");
    out->AddNode("TOWER_CALIB_HCALOUT");
    out->AddNode("CLUSTER_HCALOUT");
    out->AddNode("TOWER_SIM_CEMC");
    out->AddNode("TOWER_RAW_CEMC");
    out->AddNode("TOWER_CALIB_CEMC");
    out->AddNode("CLUSTER_CEMC");
    out->AddNode("CLUSTER_POS_COR_CEMC");
// leave the topo cluster here in case we run this during pass3
    out->AddNode("TOPOCLUSTER_ALLCALO");
    out->AddNode("TOPOCLUSTER_EMCAL");
    out->AddNode("TOPOCLUSTER_HCAL");
    out->AddNode("GlobalVertexMap");
    if (Enable::DSTOUT_COMPRESS) DstCompress(out);
    se->registerOutputManager(out);
  }
  //-----------------
  // Event processing
  //-----------------
  if (Enable::DISPLAY)
  {
    DisplayOn();

    gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
    gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");

    cout << "-------------------------------------------------" << endl;
    cout << "You are in event display mode. Run one event with" << endl;
    cout << "se->run(1)" << endl;
    cout << "Run Geant4 command with following examples" << endl;
    gROOT->ProcessLine("displaycmd()");

    return 0;
  }

  string outputfile = outputFile2 + "_g4cemc_eval.root";
  string outputfile2 = outputFile2 + "_g4hcalin_eval.root";
  string outputfile3 = outputFile2 + "_g4hcalout_eval.root";


  LiteCaloEval *eval = new LiteCaloEval("CEMCEVALUATOR", "CEMC", outputfile);
  //  eval->Verbosity(verbosity);
  eval->CaloType(LiteCaloEval::CEMC);
  se->registerSubsystem(eval);

  LiteCaloEval *eval2 = new LiteCaloEval("HCALINtowslope", "HCALIN", outputfile2);
  //  eval->Verbosity(verbosity);
  eval2->CaloType(LiteCaloEval::HCALIN);
  se->registerSubsystem(eval2);


  LiteCaloEval *eval3 = new LiteCaloEval("HOUTEVALUATOR", "HCALOUT", outputfile3);
  //  eval->Verbosity(verbosity);
  eval3->CaloType(LiteCaloEval::HCALOUT);
  se->registerSubsystem(eval3);



  // if we use a negative number of events we go back to the command line here
  if (nEvents < 0)
  {
    return 0;
  }
  // if we run the particle generator and use 0 it'll run forever
  if (nEvents == 0 && !Input::HEPMC && !Input::READHITS)
  {
    cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
    cout << "it will run forever, so I just return without running anything" << endl;
    return 0;
  }

  se->run(nEvents);

  //-----
  // Exit
  //-----

  se->End();
  se->PrintTimer();
  se->PrintMemoryTracker();
  std::cout << "All done" << std::endl;
  delete se;
  if (Enable::PRODUCTION)
  {
    Production_MoveOutput();
  }

  gSystem->Exit(0);
  return 0;
}
#endif
