//Track + Calo matching for TPC Drift Velocity Calibration
//authors: Xudong Yu <xyu3@bnl.gov>

#include <G4_ActsGeom.C>
#include <G4_Magnet.C>
#include <GlobalVariables.C>
#include <Trkr_Clustering.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>
#include <Trkr_Reco.C>
#include <G4_Global.C>

#include <caloreco/RawClusterBuilderTopo.h>

#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <tpcdvcalib/TrackToCalo.h>

#include <trackreco/AzimuthalSeeder.h>
#include <trackreco/PHActsTrackProjection.h>
#include <trackingdiagnostics/TrackResiduals.h>

#include <trackbase_historic/SvtxTrack.h>

#include <stdio.h>
#include <iostream>
#include <filesystem>

// cppcheck-suppress unknownMacro
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libTpcDVCalib.so)
R__LOAD_LIBRARY(libcalo_reco.so)

using namespace std;

namespace fs = std::filesystem;

std::string GetFirstLine(std::string listname);
bool is_directory_empty(const fs::path& dir_path);

void Fun4All_FieldOnAllTrackersCalos(
    const int nEvents = 10,
    vector<string> myInputLists = {
        "run46730_0000_trkr.txt",
        "run46730_calo.list"}, 
    std::string outDir = "./",
    bool doTpcOnlyTracking = true,
    bool doDriftVelocityCorr = true,
    int DriftVelocityCorrTag = 0,
    bool doEMcalRadiusCorr = true,
    const bool convertSeeds = false)
{
  int verbosity = 0;

  gSystem->Load("libg4dst.so");

  std::string firstfile = GetFirstLine(myInputLists[0]);
  if (*firstfile.c_str() == '\0') return;
  std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(firstfile);
  int runnumber = runseg.first;
  int segment = runseg.second;

  auto rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", 6);

  Enable::DSTOUT = false;

  TpcReadoutInit( runnumber );
  std::cout<< " run: " << runnumber
	   << " samples: " << TRACKING::reco_tpc_maxtime_sample
	   << " pre: " << TRACKING::reco_tpc_time_presample
	   << " vdrift: " << G4TPC::tpc_drift_velocity_reco
	   << std::endl;

  //Create the server
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  std::cout << "Including " << myInputLists.size() << " files." << std::endl;

  //Add all required input files
  for (unsigned int i = 0; i < myInputLists.size(); ++i)
  {
    Fun4AllInputManager *infile = new Fun4AllDstInputManager("DSTin_" + to_string(i));
    std::cout << "Including file " << myInputLists[i] << std::endl;
    //infile->AddFile(myInputLists[i]);
    infile->AddListFile(myInputLists[i]);
    se->registerInputManager(infile);
  }

  string outputAnaFileName = "TrackCalo_" + to_string(segment) + "_ana.root";
  string outputDstFileName = "TrackCalo_" + to_string(segment) + "_dst.root";

  string outputRecoDir = outDir + "/inReconstruction/" + to_string(runnumber) + "/";
  string makeDirectory = "mkdir -p " + outputRecoDir;
  system(makeDirectory.c_str());
  string outputAnaFile = outputRecoDir + outputAnaFileName;
  string outputDstFile = outputRecoDir + outputDstFileName;

  std::cout << "Reco ANA file: " << outputAnaFile << std::endl;
  std::cout << "Dst file: " << outputDstFile << std::endl;

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  double dz_separation_0;
  if (doDriftVelocityCorr)
  {
    if (DriftVelocityCorrTag==0)
    {
      // recommened drift velocity correction for Ar/CF4/N2 gas
      // tpc drift velocity = 0.00701 cm/ns
      dz_separation_0 = -0.14*2*105;
    }
    else if (DriftVelocityCorrTag==1)
    {
      // recommened drift velocity correction for Ar/CF4/ISO gas in early time, like 6x6, 28x28, 56x56 before 08/09/2024
      // tpc drift velocity = 0.00731 cm/ns
      dz_separation_0 = -0.10332944*2*105;
    }
    else if (DriftVelocityCorrTag==2)
    {
      // recommened drift velocity correction for Ar/CF4/ISO gas in 111x111 run after 08/09/2024
      // tpc drift velocity = 0.00625 cm/ns
      dz_separation_0 = -0.23335280*2*105;
    }
    else if (DriftVelocityCorrTag==3)
    {
      // recommened drift velocity correction for Ar/CF4/ISO gas in 111x111 run after 08/013/2024
      // tpc drift velocity = 0.00710 cm/ns
      dz_separation_0 = -0.12908879*2*105;
    }
    else
    {
      // for other cases which has not already been calibrated
      dz_separation_0 = 0;
    }
  }
  else
  {
    dz_separation_0 = 0;
  }
  G4TPC::tpc_drift_velocity_reco = (8.0 / 1000) * 107.0 / 105.0 * (1+ dz_separation_0 / 2. / 105.); //cm/ns
  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = true;
  //to turn on the default static corrections, enable the two lines below
  //G4TPC::ENABLE_STATIC_CORRECTIONS = true;
  //G4TPC::DISTORTIONS_USE_PHI_AS_RADIANS = false;
  //G4MAGNET::magfield = "0.01";
  G4MAGNET::magfield_tracking = G4MAGNET::magfield;
  G4MAGNET::magfield_rescale = 1;

  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;

  TRACKING::pp_mode = false;

  TrackingInit();
  if(doTpcOnlyTracking)
  {
    Tpc_HitUnpacking();
  }
  else
  {
    Mvtx_HitUnpacking();
    Intt_HitUnpacking();
    Tpc_HitUnpacking();
    Micromegas_HitUnpacking();
  }

  Mvtx_Clustering();
  Intt_Clustering();

  auto tpcclusterizer = new TpcClusterizer;
  tpcclusterizer->Verbosity(verbosity);
  tpcclusterizer->set_do_hit_association(G4TPC::DO_HIT_ASSOCIATION);
  tpcclusterizer->set_rawdata_reco();
  tpcclusterizer->set_do_sequential(true);
  se->registerSubsystem(tpcclusterizer);

  Micromegas_Clustering();

  /*
   * Begin Track Seeding
   */

  /*
   * Silicon Seeding
   */
  auto silicon_Seeding = new PHActsSiliconSeeding;
  silicon_Seeding->Verbosity(verbosity);
  silicon_Seeding->searchInIntt();
  silicon_Seeding->setinttRPhiSearchWindow(0.4);
  silicon_Seeding->setinttZSearchWindow(1.6);
  silicon_Seeding->seedAnalysis(false);
  se->registerSubsystem(silicon_Seeding);

  auto merger = new PHSiliconSeedMerger;
  merger->Verbosity(verbosity);
  se->registerSubsystem(merger);

  /*
   * Tpc Seeding
   */
  auto seeder = new PHCASeeding("PHCASeeding");
  double fieldstrength = std::numeric_limits<double>::quiet_NaN();  // set by isConstantField if constant
  bool ConstField = isConstantField(G4MAGNET::magfield_tracking, fieldstrength);
  if (ConstField)
  {
    seeder->useConstBField(true);
    seeder->constBField(fieldstrength);
  }
  else
  {
    seeder->set_field_dir(-1 * G4MAGNET::magfield_rescale);
    seeder->useConstBField(false);
    seeder->magFieldFile(G4MAGNET::magfield_tracking);  // to get charge sign right
  }
  seeder->Verbosity(verbosity);
  seeder->SetLayerRange(7, 55);
  seeder->SetSearchWindow(2.,0.05); // z-width and phi-width, default in macro at 1.5 and 0.05
  seeder->SetClusAdd_delta_window(3.0,0.06); //  (0.5, 0.005) are default; sdzdr_cutoff, d2/dr2(phi)_cutoff
  //seeder->SetNClustersPerSeedRange(4,60); // default is 6, 6
  seeder->SetMinHitsPerCluster(0);
  seeder->SetMinClustersPerTrack(3);
  seeder->useFixedClusterError(true);
  seeder->set_pp_mode(TRACKING::pp_mode);
  se->registerSubsystem(seeder);

  // expand stubs in the TPC using simple kalman filter
  auto cprop = new PHSimpleKFProp("PHSimpleKFProp");
  cprop->set_field_dir(G4MAGNET::magfield_rescale);
  if (ConstField)
  {
    cprop->useConstBField(true);
    cprop->setConstBField(fieldstrength);
  }
  else
  {
    cprop->magFieldFile(G4MAGNET::magfield_tracking);
    cprop->set_field_dir(-1 * G4MAGNET::magfield_rescale);
  }
  cprop->useFixedClusterError(true);
  cprop->set_max_window(5.);
  cprop->Verbosity(verbosity);
  cprop->set_pp_mode(TRACKING::pp_mode);
  se->registerSubsystem(cprop);

  /*
   * Track Matching between silicon and TPC
   */
  // The normal silicon association methods
  // Match the TPC track stubs from the CA seeder to silicon track stubs from PHSiliconTruthTrackSeeding
  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(verbosity);
  silicon_match->set_x_search_window(2.);
  silicon_match->set_y_search_window(2.);
  silicon_match->set_z_search_window(5.);
  silicon_match->set_phi_search_window(0.2);
  silicon_match->set_eta_search_window(0.1);
  silicon_match->set_use_old_matching(true);
  silicon_match->set_pp_mode(true);
  se->registerSubsystem(silicon_match);

  // Match TPC track stubs from CA seeder to clusters in the micromegas layers
  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(verbosity);
  mm_match->set_rphi_search_window_lyr1(0.4);
  mm_match->set_rphi_search_window_lyr2(13.0);
  mm_match->set_z_search_window_lyr1(26.0);
  mm_match->set_z_search_window_lyr2(0.4);

  mm_match->set_min_tpc_layer(38);             // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(false);  // used for tuning search windows only
  se->registerSubsystem(mm_match);

  /*
   * End Track Seeding
   */

  /*
   * Either converts seeds to tracks with a straight line/helix fit
   * or run the full Acts track kalman filter fit
   */
  if (G4TRACKING::convert_seeds_to_svtxtracks)
  {
    auto converter = new TrackSeedTrackMapConverter;
    // Option to use TpcTrackSeedContainer or SvtxTrackSeeds
    // can be set to SiliconTrackSeedContainer for silicon-only track fit
    if (doTpcOnlyTracking)
    {
      converter->setTrackSeedName("TpcTrackSeedContainer");
    }
    else
    {
      converter->setTrackSeedName("SvtxTrackSeeds");
    }
    converter->setFieldMap(G4MAGNET::magfield_tracking);
    converter->Verbosity(verbosity);
    se->registerSubsystem(converter);
  }
  else
  {
    auto deltazcorr = new PHTpcDeltaZCorrection;
    deltazcorr->Verbosity(verbosity);
    se->registerSubsystem(deltazcorr);

    // perform final track fit with ACTS
    auto actsFit = new PHActsTrkFitter;
    actsFit->Verbosity(verbosity);
    actsFit->commissioning(G4TRACKING::use_alignment);
    // in calibration mode, fit only Silicons and Micromegas hits
    if (!doTpcOnlyTracking)
    {
      actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
      actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    }
    actsFit->set_pp_mode(TRACKING::pp_mode);
    actsFit->set_use_clustermover(true);  // default is true for now
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(actsFit);
  }

  PHSimpleVertexFinder *finder = new PHSimpleVertexFinder;
  finder->Verbosity(verbosity);
  finder->setDcaCut(0.5);
  finder->setTrackPtCut(-99999.);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(1000000000);
  if (!doTpcOnlyTracking)
  {
    finder->setRequireMVTX(true);
    finder->setNmvtxRequired(3);
  }
  else
  {
    finder->setRequireMVTX(false);
  }
  finder->setOutlierPairCut(0.1);
  se->registerSubsystem(finder);

  Global_Reco();

  auto projection = new PHActsTrackProjection("CaloProjection");
  float new_cemc_rad = 100.70;//(1-(-0.077))*93.5 recommended cemc radius
  if (doEMcalRadiusCorr)
  {
    projection->setLayerRadius(SvtxTrack::CEMC, new_cemc_rad);
  }
  se->registerSubsystem(projection);

  RawClusterBuilderTopo* ClusterBuilder1 = new RawClusterBuilderTopo("EMcalRawClusterBuilderTopo1");
  ClusterBuilder1->Verbosity(verbosity);
  ClusterBuilder1->set_nodename("TOPOCLUSTER_EMCAL");
  ClusterBuilder1->set_enable_HCal(false);
  ClusterBuilder1->set_enable_EMCal(true);
  //ClusterBuilder1->set_noise(0.0025, 0.006, 0.03);
  ClusterBuilder1->set_noise(0.01, 0.03, 0.03);
  ClusterBuilder1->set_significance(4.0, 2.0, 1.0);
  ClusterBuilder1->allow_corner_neighbor(true);
  ClusterBuilder1->set_do_split(true);
  ClusterBuilder1->set_minE_local_max(1.0, 2.0, 0.5);
  ClusterBuilder1->set_R_shower(0.025);
  se->registerSubsystem(ClusterBuilder1);

  TrackToCalo *ttc = new TrackToCalo("Tracks_And_Calo", outputAnaFile);
  ttc->EMcalRadiusUser(true);
  ttc->setEMcalRadius(new_cemc_rad);
  se->registerSubsystem(ttc);

  if (Enable::DSTOUT)
  {
    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outputDstFile);
    out->StripNode("TPC");
    out->StripNode("Sync");
    out->StripNode("MBD");
    out->StripNode("ZDC");
    out->StripNode("SEPD");
    out->StripNode("HCALIN");
    out->StripNode("HCALOUT");
    out->StripNode("alignmentTransformationContainer");
    out->StripNode("alignmentTransformationContainerTransient");
    out->StripNode("SiliconTrackSeedContainer");
    out->StripNode("RUN");
    out->SaveRunNode(0);
    se->registerOutputManager(out);
  }

  se->run(nEvents);
  se->End();
  se->PrintTimer();

  ifstream file_ana(outputAnaFile.c_str(), ios::binary | ios::ate);
  if (file_ana.good() && (file_ana.tellg() > 100))
  {
    string outputRecoDirMove = outDir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputAnaFile + " " + outDir + "/Reconstructed/" + to_string(runnumber);
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  ifstream file_dst(outputDstFile.c_str(), ios::binary | ios::ate);
  if (file_dst.good() && (file_dst.tellg() > 100))
  {
    string outputRecoDirMove = outDir + "/Reconstructed/" + to_string(runnumber) + "/";
    string makeDirectoryMove = "mkdir -p " + outputRecoDirMove;
    system(makeDirectoryMove.c_str());
    string moveOutput = "mv " + outputDstFile + " " + outDir + "/Reconstructed/" + to_string(runnumber);
    std::cout << "moveOutput: " << moveOutput << std::endl;
    system(moveOutput.c_str());
  }

  delete se;
  std::cout << "All done" << std::endl;
  gSystem->Exit(0);

  return;
}

std::string GetFirstLine(std::string listname)
{
  std::ifstream file(listname);

  std::string firstLine = "";
  if (file.is_open()) {
      if (std::getline(file, firstLine)) {
          std::cout << "First Line: " << firstLine << std::endl;
      } else {
          std::cerr << "Unable to read first line of file" << std::endl;
      }
      file.close();
  } else {
      std::cerr << "Unable to open file" << std::endl;
  }
  return firstLine;
}

bool is_directory_empty(const fs::path& dir_path) {
    if (fs::exists(dir_path) && fs::is_directory(dir_path)) {
        return fs::directory_iterator(dir_path) == fs::directory_iterator();
    }
    return false;
}
