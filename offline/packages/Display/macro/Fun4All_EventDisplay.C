
int Fun4All_EventDisplay(
		       const int nEvents = 10,
		       const char * inputFile = "/sphenix/sim//sim01/production/2016-07-21/single_particle/spacal2d/fieldmap/G4Hits_sPHENIX_e-_eta0_8GeV-0002.root",
		       const char * outputFile = "G4sPHENIXCells.root",
           const char * embed_input_file = "/sphenix/sim/sim01/production/2016-07-12/sHijing/spacal2d/G4Hits_sPHENIX_sHijing-0-4.4fm.list"
		       )
{
  //===============
  // Input options
  //===============

  // Either:
  // read previously generated g4-hits files, in this case it opens a DST and skips
  // the simulations step completely. The G4Setup macro is only loaded to get information
  // about the number of layers used for the cell reco code
  //
  // In case reading production output, please double check your G4Setup_sPHENIX.C and G4_*.C consistent with those in the production macro folder
  // E.g. /sphenix/sim//sim01/production/2016-07-21/single_particle/spacal2d/
  const bool readhits = false;
  // Or:
  // read files in HepMC format (typically output from event generators like hijing or pythia)
  const bool readhepmc = false; // read HepMC files
  // Or:
  // Use particle generator
  const bool runpythia8 = true;
  const bool runpythia6 = false;
  // And
  // Further choose to embed newly simulated events to a previous simulation. Not compatible with `readhits = true`
  // In case embedding into a production output, please double check your G4Setup_sPHENIX.C and G4_*.C consistent with those in the production macro folder
  // E.g. /sphenix/sim//sim01/production/2016-07-21/single_particle/spacal2d/
  const bool do_embedding = false;

  //======================
  // What to run
  //======================

  bool do_bbc =true;
  bool do_pipe = false;
  
  bool do_svtx = true;
  bool do_svtx_cell = true;
  bool do_svtx_track = true;
  bool do_svtx_eval = true;

  bool do_preshower = false;
  
  bool do_cemc = true;
  bool do_cemc_cell = true;
  bool do_cemc_twr = true;
  bool do_cemc_cluster = true;
  bool do_cemc_eval = false;

  bool do_hcalin = true;
  bool do_hcalin_cell = true;
  bool do_hcalin_twr = true;
  bool do_hcalin_cluster = true;
  bool do_hcalin_eval = false;

  bool do_magnet = true;
  
  bool do_hcalout = true;
  bool do_hcalout_cell = true;
  bool do_hcalout_twr = true;
  bool do_hcalout_cluster = true;
  bool do_hcalout_eval = false;
  
  bool do_global = true;
  bool do_global_fastsim = false;
  
  bool do_jet_reco = true;
  bool do_jet_eval = false;

  bool do_dst_compress = false;

  //Option to convert DST to human command readable TTree for quick poke around the outputs
  bool do_DSTReader = false;
  //---------------
  // Load libraries
  //---------------

  gSystem->Load("libfun4all.so");
  gSystem->Load("libphgeom.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libphhepmc.so");
  gSystem->Load("libg4testbench.so");
  gSystem->Load("libg4hough.so");
  gSystem->Load("libcemc.so");
  gSystem->Load("libg4eval.so");

  gSystem->ListLibraries();
  // establish the geometry and reconstruction setup
  gROOT->LoadMacro("G4Setup_sPHENIX.C");
  G4Init(do_svtx,do_preshower,do_cemc,do_hcalin,do_magnet,do_hcalout,do_pipe);

  int absorberactive = 1; // set to 1 to make all absorbers active volumes
  //  const string magfield = "1.5"; // if like float -> solenoidal field in T, if string use as fieldmap name (including path)
  const string magfield = "/phenix/upgrades/decadal/fieldmaps/sPHENIX.2d.root"; // if like float -> solenoidal field in T, if string use as fieldmap name (including path)
  const float magfield_rescale = 1.4/1.5; // scale the map to a 1.4 T field

  //---------------
  // Fun4All server
  //---------------

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(2);
  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();
  // By default every random number generator uses
  // PHRandomSeed() which reads /dev/urandom to get its seed
  // if the RANDOMSEED flag is set its value is taken as seed
  // You ca neither set this to a random value using PHRandomSeed()
  // which will make all seeds identical (not sure what the point of
  // this would be:
  //  rc->set_IntFlag("RANDOMSEED",PHRandomSeed());
  // or set it to a fixed value so you can debug your code
  // rc->set_IntFlag("RANDOMSEED", 12345);

  //-----------------
  // Event generation
  //-----------------

  if (readhits)
    {
      // Get the hits from a file
      // The input manager is declared later

      if (do_embedding)
       {
         cout <<"Do not support read hits and embed background at the same time."<<endl;
         exit(1);
       }

    }
  else if (readhepmc)
    {
      // this module is needed to read the HepMC records into our G4 sims
      // but only if you read HepMC input files
      HepMCNodeReader *hr = new HepMCNodeReader();
      se->registerSubsystem(hr);
    }
  else if (runpythia8)
    {
      gSystem->Load("/afs/rhic.bnl.gov/sphenix/sys/x8664_sl6/new.2/lib/libPHPythia8.so");
      
      PHPythia8* pythia8 = new PHPythia8();
      // see coresoftware/generators/PHPythia8 for example config
      pythia8->set_config_file("phpythia8.cfg"); 
      se->registerSubsystem(pythia8);

      HepMCNodeReader *hr = new HepMCNodeReader();
      se->registerSubsystem(hr);
    }
  else if (runpythia6)
    {
      gSystem->Load("libPHPythia6.so");

      PHPythia6 *pythia6 = new PHPythia6();
      pythia6->set_config_file("phpythia6.cfg");
      se->registerSubsystem(pythia6);

      HepMCNodeReader *hr = new HepMCNodeReader();
      se->registerSubsystem(hr);
    }
  else
    {
      // toss low multiplicity dummy events
      PHG4SimpleEventGenerator *gen = new PHG4SimpleEventGenerator();
      gen->add_particles("e-",1); // mu+,e+,proton,pi+,Upsilon
      // gen->add_particles("e+",5); // mu-,e-,anti_proton,pi-
      if (readhepmc || do_embedding) {
	gen->set_reuse_existing_vertex(true);
	gen->set_existing_vertex_offset_vector(0.0,0.0,0.0);
      } else {
	gen->set_vertex_distribution_function(PHG4SimpleEventGenerator::Uniform,
					       PHG4SimpleEventGenerator::Uniform,
					       PHG4SimpleEventGenerator::Uniform);
	gen->set_vertex_distribution_mean(0.0,0.0,0.0);
	gen->set_vertex_distribution_width(0.0,0.0,5.0);
      }
      gen->set_vertex_size_function(PHG4SimpleEventGenerator::Uniform);
      gen->set_vertex_size_parameters(0.0,0.0);
      gen->set_eta_range(-0.5, 0.5);
      gen->set_phi_range(-1.0*TMath::Pi(), 1.0*TMath::Pi());
      gen->set_pt_range(0.1, 10.0);
      gen->Embed(1);
      gen->Verbosity(0);
      se->registerSubsystem(gen);
    }

  if (!readhits)
    {
      //---------------------
      // Detector description
      //---------------------

      G4Setup(absorberactive, magfield, TPythia6Decayer::kAll,
	      do_svtx, do_preshower, do_cemc, do_hcalin, do_magnet, do_hcalout, do_pipe, magfield_rescale);
    }

  //---------
  // BBC Reco
  //---------
  
  if (do_bbc) 
    {
      gROOT->LoadMacro("G4_Bbc.C");
      BbcInit();
      Bbc_Reco();
    }
  //------------------
  // Detector Division
  //------------------

  if (do_svtx_cell) Svtx_Cells();

  if (do_cemc_cell) CEMC_Cells();

  if (do_hcalin_cell) HCALInner_Cells();

  if (do_hcalout_cell) HCALOuter_Cells();

  //-----------------------------
  // CEMC towering and clustering
  //-----------------------------

  if (do_cemc_twr) CEMC_Towers();
  if (do_cemc_cluster) CEMC_Clusters();

  //-----------------------------
  // HCAL towering and clustering
  //-----------------------------
  
  if (do_hcalin_twr) HCALInner_Towers();
  if (do_hcalin_cluster) HCALInner_Clusters();

  if (do_hcalout_twr) HCALOuter_Towers();
  if (do_hcalout_cluster) HCALOuter_Clusters();

  if (do_dst_compress) ShowerCompress();

  //--------------
  // SVTX tracking
  //--------------

  if (do_svtx_track) Svtx_Reco();

  //-----------------
  // Global Vertexing
  //-----------------

  if (do_global) 
    {
      gROOT->LoadMacro("G4_Global.C");
      Global_Reco();
    }

  else if (do_global_fastsim) 
    {
      gROOT->LoadMacro("G4_Global.C");
      Global_FastSim();
    }  

  //---------
  // Jet reco
  //---------

  if (do_jet_reco) 
    {
      gROOT->LoadMacro("G4_Jets.C");
      Jet_Reco();
    }
  //----------------------
  // Simulation evaluation
  //----------------------

  if (do_svtx_eval) Svtx_Eval("g4svtx_eval.root");

  if (do_cemc_eval) CEMC_Eval("g4cemc_eval.root");

  if (do_hcalin_eval) HCALInner_Eval("g4hcalin_eval.root");

  if (do_hcalout_eval) HCALOuter_Eval("g4hcalout_eval.root");

  if (do_jet_eval) Jet_Eval("g4jet_eval.root");



  //-------------- 
  // IO management
  //--------------

  if (readhits)
    {
      // Hits file
      Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
      hitsin->fileopen(inputFile);
      se->registerInputManager(hitsin);
    }
  if (do_embedding)
    {
      if (embed_input_file == NULL)
        {
          cout << "Missing embed_input_file! Exit";
          exit(3);
        }

      Fun4AllDstInputManager *in1 = new Fun4AllNoSyncDstInputManager("DSTinEmbed");
      //      in1->AddFile(embed_input_file); // if one use a single input file
      in1->AddListFile(embed_input_file); // RecommendedL: if one use a text list of many input files
      se->registerInputManager(in1);
    }
  if (readhepmc)
    {
      Fun4AllInputManager *in = new Fun4AllHepMCInputManager( "DSTIN");
      se->registerInputManager( in );
      se->fileopen( in->Name().c_str(), inputFile );
    }
  else
    {
      // for single particle generators we just need something which drives
      // the event loop, the Dummy Input Mgr does just that
      Fun4AllInputManager *in = new Fun4AllDummyInputManager( "JADE");
      se->registerInputManager( in );
    }

  if (do_DSTReader)
    {
      //Convert DST to human command readable TTree for quick poke around the outputs
      gROOT->LoadMacro("G4_DSTReader.C");

      G4DSTreader( outputFile, //
          /*int*/ absorberactive ,
          /*bool*/ do_svtx ,
          /*bool*/ do_preshower ,
          /*bool*/ do_cemc ,
          /*bool*/ do_hcalin ,
          /*bool*/ do_magnet ,
          /*bool*/ do_hcalout ,
          /*bool*/ do_cemc_twr ,
          /*bool*/ do_hcalin_twr ,
          /*bool*/ do_magnet  ,
          /*bool*/ do_hcalout_twr
          );
    }

  // Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outputFile);
  // if (do_dst_compress) DstCompress(out);
  // se->registerOutputManager(out);


  //-----------------
  // Event processing
  //-----------------
  if (nEvents < 0)
    {
      return;
    }
  // if we run the particle generator and use 0 it'll run forever
  if (nEvents == 0 && !readhits && !readhepmc)
    {
      cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
      cout << "it will run forever, so I just return without running anything" << endl;
      return;
    }

    //---------------
    // Event display
    //---------------

     gSystem->Exec("/bin/env");
     gSystem->Load("libpheve.so");
     gSystem->Load("libPHBFieldMap.so");
     //          	    width height use_fieldmap use_geofile field-map name    geo-fiel name
     PHEventDisplay* event_display
        = new PHEventDisplay(1920,860,     true,      false,     "sPHENIX.2d.root", "sphenix_maps+tpc_geo.root");
     // location of geometry file : /gpfs/mnt/gpfs02/sphenix/user/shlee/svtx/stage1_jobs
     event_display->set_jet_pt_threshold(10.);//GeV/c
     event_display->set_jet_e_scale(30.);//GeV
     event_display->set_calo_e_threshold(0.0009);//GeV
     event_display->set_svtx_on(true);
     event_display->set_cemc_on(true);
     event_display->set_hcalin_on(false);
     event_display->set_hcalout_on(false);
     event_display->set_jet_on(false);
     event_display->set_truth_on(false);
     event_display->set_verbosity(2);

     se->registerSubsystem(event_display);
     //event_display->start_rotation();
     timer = new TTimer();
     timer->Connect("Timeout()", "PHEventDisplay", event_display,"run_evt_in_thread()");
     timer->Start(100000,kFALSE);

}
