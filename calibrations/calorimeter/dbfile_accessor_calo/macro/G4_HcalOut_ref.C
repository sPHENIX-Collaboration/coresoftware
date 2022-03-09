#ifndef MACRO_G4HCALOUTREF_C
#define MACRO_G4HCALOUTREF_C

#include <GlobalVariables.C>
#include <QA.C>

#include <g4calo/HcalRawTowerBuilder.h>
#include <g4calo/RawTowerDigitizer.h>

#include <g4detectors/PHG4HcalCellReco.h>
#include <g4detectors/PHG4OuterHcalSubsystem.h>

#include <g4eval/CaloEvaluator.h>

#include <g4main/PHG4Reco.h>

#include <caloreco/RawClusterBuilderGraph.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawTowerCalibration.h>
#include <qa_modules/QAG4SimulationCalorimeter.h>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libqa_modules.so)

namespace Enable
{
  bool HCALOUT = false;
  bool HCALOUT_ABSORBER = false;
  bool HCALOUT_OVERLAPCHECK = false;
  bool HCALOUT_CELL = false;
  bool HCALOUT_TOWER = false;
  bool HCALOUT_CLUSTER = false;
  bool HCALOUT_EVAL = false;
  bool HCALOUT_QA = false;
  int HCALOUT_VERBOSITY = 0;
}  // namespace Enable

namespace G4HCALOUT
{
  double outer_radius = 264.71;
  double size_z = 304.91 * 2;

  // Digitization (default photon digi):
  RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kSimple_photon_digitization;
  // directly pass the energy of sim tower to digitized tower
  // kNo_digitization
  // simple digitization with photon statistics, single amplitude ADC conversion and pedestal
  // kSimple_photon_digitization
  // digitization with photon statistics on SiPM with an effective pixel N, ADC conversion and pedestal
  // kSiPM_photon_digitization

  enum enu_HCalOut_clusterizer
  {
    kHCalOutGraphClusterizer,
    kHCalOutTemplateClusterizer
  };

  //! template clusterizer, RawClusterBuilderTemplate, as developed by Sasha Bazilevsky
  enu_HCalOut_clusterizer HCalOut_clusterizer = kHCalOutTemplateClusterizer;
  //! graph clusterizer, RawClusterBuilderGraph
  //enu_HCalOut_clusterizer HCalOut_clusterizer = kHCalOutGraphClusterizer;
}  // namespace G4HCALOUT

// Init is called by G4Setup.C
void HCalOuterInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4HCALOUT::outer_radius);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, G4HCALOUT::size_z / 2.);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -G4HCALOUT::size_z / 2.);
}

double HCalOuter(PHG4Reco *g4Reco,
                 double radius,
                 const int crossings)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::HCALOUT_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::HCALOUT_OVERLAPCHECK;
  int verbosity = std::max(Enable::VERBOSITY, Enable::HCALOUT_VERBOSITY);

  PHG4OuterHcalSubsystem *hcal = new PHG4OuterHcalSubsystem("HCALOUT");
  // hcal->set_double_param("inner_radius", 183.3);
  //-----------------------------------------
  // the light correction can be set in a single call
  // hcal->set_double_param("light_balance_inner_corr", NAN);
  // hcal->set_double_param("light_balance_inner_radius", NAN);
  // hcal->set_double_param("light_balance_outer_corr", NAN);
  // hcal->set_double_param("light_balance_outer_radius", NAN);
  // hcal->set_double_param("magnet_cutout_radius", 195.31);
  // hcal->set_double_param("magnet_cutout_scinti_radius", 195.96);
  // hcal->SetLightCorrection(NAN,NAN,NAN,NAN);
  //-----------------------------------------
  // hcal->set_double_param("outer_radius", G4HCALOUT::outer_radius);
  // hcal->set_double_param("place_x", 0.);
  // hcal->set_double_param("place_y", 0.);
  // hcal->set_double_param("place_z", 0.);
  // hcal->set_double_param("rot_x", 0.);
  // hcal->set_double_param("rot_y", 0.);
  // hcal->set_double_param("rot_z", 0.);
  // hcal->set_double_param("scinti_eta_coverage", 1.1);
  // hcal->set_double_param("scinti_gap", 0.85);
  // hcal->set_double_param("scinti_gap_neighbor", 0.1);
  // hcal->set_double_param("scinti_inner_radius",183.89);
  // hcal->set_double_param("scinti_outer_radius",263.27);
  // hcal->set_double_param("scinti_tile_thickness", 0.7);
  // hcal->set_double_param("size_z", G4HCALOUT::size_z);
  // hcal->set_double_param("steplimits", NAN);
  // hcal->set_double_param("tilt_angle", -11.23);

  // hcal->set_int_param("light_scint_model", 1);
  // hcal->set_int_param("magnet_cutout_first_scinti", 8);
  // hcal->set_int_param("ncross", 0);
  // hcal->set_int_param("n_towers", 64);
  // hcal->set_int_param("n_scinti_plates_per_tower", 5);
  // hcal->set_int_param("n_scinti_tiles", 12);

  // hcal->set_string_param("material", "Steel_1006");

  hcal->SetActive();
  hcal->SuperDetector("HCALOUT");
  if (AbsorberActive)
  {
    hcal->SetAbsorberActive();
  }
  hcal->OverlapCheck(OverlapCheck);
  g4Reco->registerSubsystem(hcal);

  radius = hcal->get_double_param("outer_radius");

  radius += no_overlapp;

  return radius;
}

void HCALOuter_Cells()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::HCALOUT_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  PHG4HcalCellReco *hc = new PHG4HcalCellReco("HCALOUT_CELLRECO");
  hc->Detector("HCALOUT");
  //  hc->Verbosity(2);
  // check for energy conservation - needs modified "infinite" timing cuts
  // 0-999999999
  //  hc->checkenergy();
  // timing cuts with their default settings
  // hc->set_double_param("tmin",0.);
  // hc->set_double_param("tmax",60.0);
  // or all at once:
  // hc->set_timing_window(0.0,60.0);
  // this sets all cells to a fixed energy for debugging
  // hc->set_fixed_energy(1.);
  se->registerSubsystem(hc);

  return;
}

void HCALOuter_Towers()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::HCALOUT_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  HcalRawTowerBuilder *TowerBuilder = new HcalRawTowerBuilder("HcalOutRawTowerBuilder");
  TowerBuilder->Detector("HCALOUT");
  TowerBuilder->set_sim_tower_node_prefix("SIM");
  // this sets specific decalibration factors
  // for a given cell
  // TowerBuilder->set_cell_decal_factor(1,10,0.1);
  // for a whole tower
  // TowerBuilder->set_tower_decal_factor(0,10,0.2);
  // TowerBuilder->set_cell_decal_factor(1,10,0.1);
  // TowerBuilder->set_tower_decal_factor(0,10,0.2);
  TowerBuilder->Verbosity(verbosity);
  se->registerSubsystem(TowerBuilder);

  // From 2016 Test beam sim
  RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("HcalOutRawTowerDigitizer");
  TowerDigitizer->Detector("HCALOUT");
  //  TowerDigitizer->set_raw_tower_node_prefix("RAW_LG");
  TowerDigitizer->set_digi_algorithm(G4HCALOUT::TowerDigi);
  TowerDigitizer->set_pedstal_central_ADC(0);
  TowerDigitizer->set_pedstal_width_ADC(1);  // From Jin's guess. No EMCal High Gain data yet! TODO: update
  TowerDigitizer->set_photonelec_ADC(16. / 5.);
  TowerDigitizer->set_photonelec_yield_visible_GeV(16. / 5 / (0.2e-3));
  TowerDigitizer->set_zero_suppression_ADC(-0);  // no-zero suppression
  se->registerSubsystem(TowerDigitizer);

  const double visible_sample_fraction_HCALOUT = 3.38021e-02;  // /gpfs/mnt/gpfs04/sphenix/user/jinhuang/prod_analysis/hadron_shower_res_nightly/./G4Hits_sPHENIX_pi-_eta0_16GeV.root_qa.rootQA_Draw_HCALOUT_G4Hit.pdf

  RawTowerCalibration *TowerCalibration = new RawTowerCalibration("HcalOutRawTowerCalibration");
  TowerCalibration->Detector("HCALOUT");
  //  TowerCalibration->set_raw_tower_node_prefix("RAW_LG");
  //  TowerCalibration->set_calib_tower_node_prefix("CALIB_LG");

  //TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);

  TowerCalibration->set_calib_algorithm(RawTowerCalibration::kDbfile_tbt_gain_corr);
  TowerCalibration->set_UseConditionsDB(false);
  //  TowerCalibration->set_CalibrationFileName("HCALOUT_GainsCalib1.22.txt"); // currently to be picked up in running directory
  TowerCalibration->set_CalibrationFileName("HCALOUT_GainsCalib.txt"); // in calibrations repos, $CALIBRATIONROOT/HCALOUT/
  
  
  if (G4HCALOUT::TowerDigi == RawTowerDigitizer::kNo_digitization)
  {
    // 0.033 extracted from electron sims (edep(scintillator)/edep(total))
    TowerCalibration->set_calib_const_GeV_ADC(1. / 0.033);
  }
  else
  {
    TowerCalibration->set_calib_const_GeV_ADC(0.2e-3 / visible_sample_fraction_HCALOUT);
  }
  TowerCalibration->set_pedstal_ADC(0);
  se->registerSubsystem(TowerCalibration);

  return;
}

void HCALOuter_Clusters()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::HCALOUT_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  if (G4HCALOUT::HCalOut_clusterizer == G4HCALOUT::kHCalOutTemplateClusterizer)
  {
    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("HcalOutRawClusterBuilderTemplate");
    ClusterBuilder->Detector("HCALOUT");
    ClusterBuilder->SetCylindricalGeometry();  // has to be called after Detector()
    ClusterBuilder->Verbosity(verbosity);
    se->registerSubsystem(ClusterBuilder);
  }
  else if (G4HCALOUT::HCalOut_clusterizer == G4HCALOUT::kHCalOutGraphClusterizer)
  {
    RawClusterBuilderGraph *ClusterBuilder = new RawClusterBuilderGraph("HcalOutRawClusterBuilderGraph");
    ClusterBuilder->Detector("HCALOUT");
    ClusterBuilder->Verbosity(verbosity);
    se->registerSubsystem(ClusterBuilder);
  }
  else
  {
    cout << "HCALOuter_Clusters - unknown clusterizer setting!" << endl;
    exit(1);
  }

  return;
}

void HCALOuter_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::HCALOUT_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloEvaluator *eval = new CaloEvaluator("HCALOUTEVALUATOR", "HCALOUT", outputfile);
  eval->Verbosity(verbosity);
  se->registerSubsystem(eval);

  return;
}

void HCALOuter_QA()
{
  int verbosity = std::max(Enable::QA_VERBOSITY, Enable::HCALOUT_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();
  QAG4SimulationCalorimeter *qa = new QAG4SimulationCalorimeter("HCALOUT");
  qa->Verbosity(verbosity);
  se->registerSubsystem(qa);

  return;
}

#endif
