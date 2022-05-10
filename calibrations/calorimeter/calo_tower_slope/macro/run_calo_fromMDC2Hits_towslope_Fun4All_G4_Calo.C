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

#include <litecaloeval/LiteCaloEval.h>
#include <calib_emc_pi0/CaloCalibEmc_Pi0.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcaloCalibDBFile.so)
R__LOAD_LIBRARY(libcalibCaloEmc_pi0.so)
R__LOAD_LIBRARY(libLiteCaloEvalTowSlope.so)


// For HepMC Hijing
// try inputFile = /sphenix/sim/sim01/sphnxpro/sHijing_HepMC/sHijing_0-12fm.dat

int run_calo_fromMDC2Hits_towslope_Fun4All_G4_Calo(
    const int nEvents = 2,
    //    const string &inputFile0 = "DST_CALO_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000004-10000.root",
    //const string &inputFile1 = "DST_VERTEX_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000004-10000.root",
    const int mdc2_4_file_num = 1,
    const string &outputFile = "newoutput1_calo1",
    //      const string &inputFile0 = "/gpfs/mnt/gpfs02/sphenix/sim/sim01/sphnxpro/MDC2/sHijing_HepMC/fm_0_20/pileup/DST_CALO_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000003-00001.root",  //"DST_CALO_G4HIT_sHijing_0_12fm-0000000001-00000.root",
    //      const string &inputFile1 = "/gpfs/mnt/gpfs02/sphenix/sim/sim01/sphnxpro/MDC2/sHijing_HepMC/fm_0_20/pileup/DST_VERTEX_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000003-00001.root", //"DST_VERTEX_sHijing_0_12fm-0000000001-00000.root",
    //      const string &outputFile = "G4sPHENIX_calo1",
      const string &outdir = ".")
{
   Fun4AllServer *se = Fun4AllServer::instance();
   se->Verbosity(0);


   //   string inputFile0 = "DST_CALO_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000004-10000.root";
   //   string inputFile1 = "DST_VERTEX_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000004-10000.root";

   //There are 40,000 mdc2 production 4 files.   
   // this piece of code allows you to simply specify which of the 40k files 
   // by the single integer input (from 0-40k) mdc2_4_file_num

   string inputFile0 = "DST_CALO_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000004-";
   string inputFile1 = "DST_VERTEX_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000004-";

   
   int ynum_int = 100000+ mdc2_4_file_num;
   TString yn_tstr = "";
   yn_tstr += ynum_int;
   yn_tstr.Remove(0,1);
   inputFile0 += yn_tstr.Data();
   inputFile1 += yn_tstr.Data();

   inputFile0 += ".root";
   inputFile1 += ".root";
   
   cout << "running over these files" << endl;
   cout << inputFile0 << endl;
   cout << inputFile1 << endl;

   //const string &outputFile = "newoutput1_calo1",


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

  Enable::DSTOUT = false;
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

  string outputfile = outputFile2 + "_g4hcin_eval.root";
  string outputfile2 = outputFile2 + "_hout.root";
  string outputfile3 = outputFile2 + "_hin_eval_mod.root";
  string outputfile4 = outputFile2 + "_hout_mod.root";
  string outputfile5 = outputFile2 + "_cemc_eval.root";
  string outputfile6 = outputFile2 + "_cemceval_mod.root";


  LiteCaloEval *eval = new LiteCaloEval("HOUTEVALUATOR0", "HCALOUT", outputfile2.c_str());
  //"CEMCEVALUATOR2", "CEMC", outputfile.c_str());
  //  eval->Verbosity(verbosity);
  eval->CaloType(LiteCaloEval::HCALOUT);
  se->registerSubsystem(eval);


  // LiteCaloEval *eval3a = new LiteCaloEval("HOUTEVALUATOR", "HCALOUT", outputfile4.c_str());
  // //  eval->Verbosity(verbosity);
  // eval3a->CaloType(LiteCaloEval::HCALOUT);  //  eval->Verbosity(verbosity);
  // eval3a->set_mode(1);
  // se->registerSubsystem(eval3a);


  LiteCaloEval *eval2 = new LiteCaloEval("HINEVALUATOR2", "HCALIN", outputfile.c_str());
    //  eval->Verbosity(verbosity);
  eval2->CaloType(LiteCaloEval::HCALIN);
  se->registerSubsystem(eval2);


  // LiteCaloEval *eval4a = new LiteCaloEval("HINEVALUATOR4", "HCALIN", outputfile3.c_str());
  // //  eval->Verbosity(verbosity);
  // eval4a->CaloType(LiteCaloEval::HCALIN);  //  eval->Verbosity(verbosity);
  // eval4a->set_mode(1);
  // se->registerSubsystem(eval4a);



  LiteCaloEval *eval5 = new LiteCaloEval("CEMCEVALUATOR2", "CEMC", outputfile5.c_str());
  //  eval->Verbosity(verbosity);
  eval5->CaloType(LiteCaloEval::CEMC);
  se->registerSubsystem(eval5);


  // LiteCaloEval *eval6 = new LiteCaloEval("CEMCEVALUATOR", "CEMC", outputfile6.c_str());
  // //  eval->Verbosity(verbosity);
  // eval6->CaloType(LiteCaloEval::CEMC);
  // eval6->set_mode(1);
  // se->registerSubsystem(eval6);



  /*
  LiteCaloEval *eval = new LiteCaloEval("CEMCEVALUATOR2", "CEMC", outputfile.c_str());
  //  eval->Verbosity(verbosity);
  eval->CaloType(LiteCaloEval::CEMC);
  se->registerSubsystem(eval);

  */
  /*
  CaloCalibEmc_Pi0 *eval_pi1 = new CaloCalibEmc_Pi0("CEMC_CALIB_PI0", outputfile2);
  //  eval_pi1->set_mode(1);
  //  eval->Verbosity(verbosity);
  se->registerSubsystem(eval_pi1);
  */

  /*
  CaloCalibEmc_Pi0 *eval_pi2 = new CaloCalibEmc_Pi0("CEMC_CALIB_PI02", outputfile4);
  eval_pi2->set_mode(1);
  //  eval->Verbosity(verbosity);
  se->registerSubsystem(eval_pi2);
  */

  /*
  LiteCaloEval *eval3 = new LiteCaloEval("HOUTEVALUATOR", "HCALOUT", outputfile3.c_str());
  //  eval->Verbosity(verbosity);
  eval3->CaloType(LiteCaloEval::HCALOUT);
  //  se->registerSubsystem(eval3);
  */

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
  //  se->PrintTimer();
  //se->PrintMemoryTracker();
  std::cout << "All done" << std::endl;
  gSystem->Exit(0);
  delete se;
  /*
  if (Enable::PRODUCTION)
  {
    Production_MoveOutput();
  }
  */

  return 0;
}
#endif
