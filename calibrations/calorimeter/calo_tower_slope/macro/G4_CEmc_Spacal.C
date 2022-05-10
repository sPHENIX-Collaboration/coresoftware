#ifndef MACRO_G4CEMCSPACAL_C
#define MACRO_G4CEMCSPACAL_C

#include <GlobalVariables.C>
#include <QA.C>

#include <g4detectors/PHG4CylinderCellReco.h>
#include <g4detectors/PHG4CylinderGeom_Spacalv1.h>
#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4detectors/PHG4FullProjSpacalCellReco.h>
#include <g4detectors/PHG4SpacalSubsystem.h>

#include <g4calo/RawTowerBuilder.h>
#include <g4calo/RawTowerDigitizer.h>

#include <g4eval/CaloEvaluator.h>

#include <g4main/PHG4Reco.h>
#include <g4main/PHG4Utils.h>

#include <caloreco/RawClusterBuilderGraph.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include <caloreco/RawTowerCalibration.h>
#include <qa_modules/QAG4SimulationCalorimeter.h>

#include <fun4all/Fun4AllServer.h>

double
CEmc_1DProjectiveSpacal(PHG4Reco *g4Reco, double radius, const int crossings);

double
CEmc_2DProjectiveSpacal(PHG4Reco *g4Reco, double radius, const int crossings);

R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4calo.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libqa_modules.so)

namespace Enable
{
  bool CEMC = false;
  bool CEMC_ABSORBER = false;
  bool CEMC_OVERLAPCHECK = false;
  bool CEMC_CELL = false;
  bool CEMC_TOWER = false;
  bool CEMC_CLUSTER = false;
  bool CEMC_EVAL = false;
  bool CEMC_QA = false;
  int CEMC_VERBOSITY = 0;
}  // namespace Enable

namespace G4CEMC
{
  int Min_cemc_layer = 1;
  int Max_cemc_layer = 1;

  // Digitization (default photon digi):
  RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kSimple_photon_digitization;
  //  RawTowerDigitizer::enu_digi_algorithm TowerDigi = RawTowerDigitizer::kNo_digitization;
  // directly pass the energy of sim tower to digitized tower
  // kNo_digitization
  // simple digitization with photon statistics, single amplitude ADC conversion and pedestal
  // kSimple_photon_digitization
  // digitization with photon statistics on SiPM with an effective pixel N, ADC conversion and pedestal
  // kSiPM_photon_digitization

  // set a default value for SPACAL configuration
  //  // 1D azimuthal projective SPACAL (fast)
  //int Cemc_spacal_configuration = PHG4CylinderGeom_Spacalv1::k1DProjectiveSpacal;
  //   2D azimuthal projective SPACAL (slow)
  int Cemc_spacal_configuration = PHG4CylinderGeom_Spacalv1::k2DProjectiveSpacal;

  enum enu_Cemc_clusterizer
  {
    kCemcGraphClusterizer,

    kCemcTemplateClusterizer
  };

  //! template clusterizer, RawClusterBuilderTemplate, as developed by Sasha Bazilevsky
  enu_Cemc_clusterizer Cemc_clusterizer = kCemcTemplateClusterizer;
  //! graph clusterizer, RawClusterBuilderGraph
  //enu_Cemc_clusterizer Cemc_clusterizer = kCemcGraphClusterizer;

}  // namespace G4CEMC

// black hole parameters are set in CEmc function
// needs a dummy argument to play with current G4Setup_sPHENIX.C
void CEmcInit(const int i = 0)
{
}

//! EMCal main setup macro
double
CEmc(PHG4Reco *g4Reco, double radius, const int crossings)
{
  if (G4CEMC::Cemc_spacal_configuration == PHG4CylinderGeom_Spacalv1::k1DProjectiveSpacal)
  {
    return CEmc_1DProjectiveSpacal(/*PHG4Reco**/ g4Reco, /*double*/ radius, /*const int */ crossings);
  }
  else if (G4CEMC::Cemc_spacal_configuration == PHG4CylinderGeom_Spacalv1::k2DProjectiveSpacal)
  {
    return CEmc_2DProjectiveSpacal(/*PHG4Reco**/ g4Reco, /*double*/ radius, /*const int */ crossings);
  }
  else
  {
    std::cout
        << "G4_CEmc_Spacal.C::CEmc - Fatal Error - unrecognized SPACAL configuration #"
        << G4CEMC::Cemc_spacal_configuration << ". Force exiting..." << std::endl;
    exit(-1);
    return 0;
  }
}

//! EMCal setup macro - 1D azimuthal projective SPACAL
double
CEmc_1DProjectiveSpacal(PHG4Reco *g4Reco, double radius, const int crossings)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::CEMC_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::CEMC_OVERLAPCHECK;

  double emc_inner_radius = 95.;  // emc inner radius from engineering drawing
  double cemcthickness = 12.7;
  double emc_outer_radius = emc_inner_radius + cemcthickness;  // outer radius

  if (radius > emc_inner_radius)
  {
    cout << "inconsistency: pstof outer radius: " << radius
         << " larger than emc inner radius: " << emc_inner_radius
         << endl;
    gSystem->Exit(-1);
  }

  //  boundary check
  if (radius > emc_inner_radius - 1.5 - no_overlapp)
  {
    cout << "G4_CEmc_Spacal.C::CEmc() - expect radius < " << emc_inner_radius - 1.5 - no_overlapp << " to install SPACAL" << endl;
    exit(1);
  }
  radius = emc_inner_radius - 1.5 - no_overlapp;

  // 1.5cm thick teflon as an approximation for EMCAl light collection + electronics (10% X0 total estimated)
  PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem("CEMC_ELECTRONICS", 0);
  cyl->SuperDetector("CEMC_ELECTRONICS");
  cyl->set_double_param("radius", radius);
  cyl->set_string_param("material", "G4_TEFLON");
  cyl->set_double_param("thickness", 1.5);
  if (AbsorberActive) cyl->SetActive();
  g4Reco->registerSubsystem(cyl);

  radius += 1.5;
  radius += no_overlapp;

  int ilayer = G4CEMC::Min_cemc_layer;
  PHG4SpacalSubsystem *cemc = new PHG4SpacalSubsystem("CEMC", ilayer);
  cemc->set_double_param("radius", emc_inner_radius);
  cemc->set_double_param("thickness", cemcthickness);

  cemc->SetActive();
  cemc->SuperDetector("CEMC");
  if (AbsorberActive) cemc->SetAbsorberActive();
  cemc->OverlapCheck(OverlapCheck);

  g4Reco->registerSubsystem(cemc);

  if (ilayer > G4CEMC::Max_cemc_layer)
  {
    cout << "layer discrepancy, current layer " << ilayer
         << " max cemc layer: " << G4CEMC::Max_cemc_layer << endl;
  }

  radius += cemcthickness;
  radius += no_overlapp;

  // 0.5cm thick Stainless Steel as an approximation for EMCAl support system
  cyl = new PHG4CylinderSubsystem("CEMC_SPT", 0);
  cyl->SuperDetector("CEMC_SPT");
  cyl->set_double_param("radius", radius);
  cyl->set_string_param("material", "SS310");  // SS310 Stainless Steel
  cyl->set_double_param("thickness", 0.5);
  if (AbsorberActive) cyl->SetActive();
  g4Reco->registerSubsystem(cyl);
  radius += 0.5;
  // this is the z extend and outer radius of the support structure and therefore the z extend
  // and radius of the surrounding black holes
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 149.47);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -149.47);
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, radius);
  radius += no_overlapp;

  return radius;
}

//! 2D full projective SPACAL
double
CEmc_2DProjectiveSpacal(PHG4Reco *g4Reco, double radius, const int crossings)
{
  bool AbsorberActive = Enable::ABSORBER || Enable::CEMC_ABSORBER;
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::CEMC_OVERLAPCHECK;

  double emc_inner_radius = 92;  // emc inner radius from engineering drawing
  double cemcthickness = 24.00000 - no_overlapp;

  //max radius is 116 cm;
  double emc_outer_radius = emc_inner_radius + cemcthickness;  // outer radius
  assert(emc_outer_radius < 116);

  if (radius > emc_inner_radius)
  {
    cout << "inconsistency: preshower radius+thickness: " << radius
         << " larger than emc inner radius: " << emc_inner_radius << endl;
    gSystem->Exit(-1);
  }

  // the radii are only to determined the thickness of the cemc
  radius = emc_inner_radius;

  // 1.5cm thick teflon as an approximation for EMCAl light collection + electronics (10% X0 total estimated)
  PHG4CylinderSubsystem *cyl = new PHG4CylinderSubsystem("CEMC_ELECTRONICS", 0);
  cyl->set_double_param("radius", radius);
  cyl->set_string_param("material", "G4_TEFLON");
  cyl->set_double_param("thickness", 1.5 - no_overlapp);
  cyl->SuperDetector("CEMC_ELECTRONICS");
  cyl->OverlapCheck(OverlapCheck);
  if (AbsorberActive) cyl->SetActive();
  g4Reco->registerSubsystem(cyl);

  radius += 1.5;
  cemcthickness -= 1.5 + no_overlapp;

  // 0.5cm thick Stainless Steel as an approximation for EMCAl support system
  cyl = new PHG4CylinderSubsystem("CEMC_SPT", 0);
  cyl->SuperDetector("CEMC_SPT");
  cyl->set_double_param("radius", radius + cemcthickness - 0.5);
  cyl->set_string_param("material", "SS310");  // SS310 Stainless Steel
  cyl->set_double_param("thickness", 0.5 - no_overlapp);
  cyl->OverlapCheck(OverlapCheck);
  if (AbsorberActive) cyl->SetActive();
  g4Reco->registerSubsystem(cyl);

  // this is the z extend and outer radius of the support structure and therefore the z extend
  // and radius of the surrounding black holes
  double sptlen = PHG4Utils::GetLengthForRapidityCoverage(radius + cemcthickness);
  BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, sptlen);
  BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -sptlen);
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, radius + cemcthickness);

  cemcthickness -= 0.5 + no_overlapp;

  int ilayer = 0;
  PHG4SpacalSubsystem *cemc;

  cemc = new PHG4SpacalSubsystem("CEMC", ilayer);

  cemc->set_int_param("virualize_fiber", 0);
  cemc->set_int_param("azimuthal_seg_visible", 1);
  cemc->set_int_param("construction_verbose", 0);
  cemc->Verbosity(0);

  cemc->UseCalibFiles(PHG4DetectorSubsystem::xml);
  cemc->SetCalibrationFileDir(string(getenv("CALIBRATIONROOT")) + string("/CEMC/Geometry_2018ProjTilted/"));
  cemc->set_double_param("radius", radius);            // overwrite minimal radius
  cemc->set_double_param("thickness", cemcthickness);  // overwrite thickness

  cemc->SetActive();
  cemc->SuperDetector("CEMC");
  if (AbsorberActive) cemc->SetAbsorberActive();
  cemc->OverlapCheck(OverlapCheck);

  g4Reco->registerSubsystem(cemc);

  if (ilayer > G4CEMC::Max_cemc_layer)
  {
    cout << "layer discrepancy, current layer " << ilayer
         << " max cemc layer: " << G4CEMC::Max_cemc_layer << endl;
  }

  radius += cemcthickness;
  radius += no_overlapp;

  return radius;
}

void CEMC_Cells()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::CEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  if (G4CEMC::Cemc_spacal_configuration == PHG4CylinderGeom_Spacalv1::k1DProjectiveSpacal)
  {
    PHG4CylinderCellReco *cemc_cells = new PHG4CylinderCellReco("CEMCCYLCELLRECO");
    cemc_cells->Detector("CEMC");
    cemc_cells->Verbosity(verbosity);
    for (int i = G4CEMC::Min_cemc_layer; i <= G4CEMC::Max_cemc_layer; i++)
    {
      //          cemc_cells->etaphisize(i, 0.024, 0.024);
      const double radius = 95;
      cemc_cells->cellsize(i, 2 * M_PI / 256. * radius, 2 * M_PI / 256. * radius);
    }
    se->registerSubsystem(cemc_cells);
  }
  else if (G4CEMC::Cemc_spacal_configuration == PHG4CylinderGeom_Spacalv1::k2DProjectiveSpacal)
  {
    PHG4FullProjSpacalCellReco *cemc_cells = new PHG4FullProjSpacalCellReco("CEMCCYLCELLRECO");
    cemc_cells->Detector("CEMC");
    cemc_cells->Verbosity(verbosity);
    cemc_cells->get_light_collection_model().load_data_file(
        string(getenv("CALIBRATIONROOT")) + string("/CEMC/LightCollection/Prototype3Module.xml"),
        "data_grid_light_guide_efficiency", "data_grid_fiber_trans");
    se->registerSubsystem(cemc_cells);
  }
  else
  {
    cout << "G4_CEmc_Spacal.C::CEmc - Fatal Error - unrecognized SPACAL configuration #"
         << G4CEMC::Cemc_spacal_configuration << ". Force exiting..." << endl;
    gSystem->Exit(-1);
    return;
  }

  return;
}

void CEMC_Towers()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::CEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  RawTowerBuilder *TowerBuilder = new RawTowerBuilder("EmcRawTowerBuilder");
  TowerBuilder->Detector("CEMC");
  TowerBuilder->set_sim_tower_node_prefix("SIM");
  TowerBuilder->Verbosity(verbosity);
  se->registerSubsystem(TowerBuilder);

  double sampling_fraction = 1;
  if (G4CEMC::Cemc_spacal_configuration == PHG4CylinderGeom_Spacalv1::k1DProjectiveSpacal)
  {
    sampling_fraction = 0.0234335;  //from production:/gpfs02/phenix/prod/sPHENIX/preCDR/pro.1-beta.3/single_particle/spacal1d/zerofield/G4Hits_sPHENIX_e-_eta0_8GeV.root
  }
  else if (G4CEMC::Cemc_spacal_configuration == PHG4CylinderGeom_Spacalv1::k2DProjectiveSpacal)
  {
    //      sampling_fraction = 0.02244; //from production: /gpfs02/phenix/prod/sPHENIX/preCDR/pro.1-beta.3/single_particle/spacal2d/zerofield/G4Hits_sPHENIX_e-_eta0_8GeV.root
    //    sampling_fraction = 2.36081e-02;  //from production: /gpfs02/phenix/prod/sPHENIX/preCDR/pro.1-beta.5/single_particle/spacal2d/zerofield/G4Hits_sPHENIX_e-_eta0_8GeV.root
    //    sampling_fraction = 1.90951e-02; // 2017 Tilt porjective SPACAL, 8 GeV photon, eta = 0.3 - 0.4
    sampling_fraction = 2e-02;  // 2017 Tilt porjective SPACAL, tower-by-tower calibration
  }
  else
  {
    std::cout
        << "G4_CEmc_Spacal.C::CEMC_Towers - Fatal Error - unrecognized SPACAL configuration #"
        << G4CEMC::Cemc_spacal_configuration << ". Force exiting..." << std::endl;
    exit(-1);
    return;
  }

  const double photoelectron_per_GeV = 500;  //500 photon per total GeV deposition

  bool doSimple = true;

  RawTowerDigitizer *TowerDigitizer = new RawTowerDigitizer("EmcRawTowerDigitizer");
  TowerDigitizer->Detector("CEMC");
  TowerDigitizer->Verbosity(verbosity);
  TowerDigitizer->set_digi_algorithm(G4CEMC::TowerDigi); 
  TowerDigitizer->set_variable_pedestal(true);  //read ped central and width from calibrations file comment next 2 lines if true
                                                //  TowerDigitizer->set_pedstal_central_ADC(0);
                                                //  TowerDigitizer->set_pedstal_width_ADC(8);  // eRD1 test beam setting
  TowerDigitizer->set_photonelec_ADC(1);        //not simulating ADC discretization error
  TowerDigitizer->set_photonelec_yield_visible_GeV(photoelectron_per_GeV / sampling_fraction);
  TowerDigitizer->set_variable_zero_suppression(true);  //read zs values from calibrations file comment next line if true
                                                        //  TowerDigitizer->set_zero_suppression_ADC(16);  // eRD1 test beam setting
  TowerDigitizer->GetParameters().ReadFromFile("CEMC", "xml", 0, 0,
                                               string(getenv("CALIBRATIONROOT")) + string("/CEMC/TowerCalibCombinedParams_2020/"));  // calibration database


  // tower digitizer settings for doing decalibration -JEF May '22

  TowerDigitizer->set_UseConditionsDB(false);

  //---------------
  // if conditions db is enabled in the line above, how to handle filename
  // needs decided see below examples  (TBD by Chris P)
  //------------------
  // Some standard decal files specification where full decal is happening
  // uses the same db file accessor api as for TowerCalib and thus format 
  // (format can change along with accessor internals, but user
  // needs to know/give a readable format file).
  // Decal tower by tow. factors in file are apply as a multiplicative 
  // factor to raw energy/adc see below about adc level
  //---
  //  TowerDigitizer->set_DoTowerDecal(true,"emcal_corr_sec12bands.root",false);
  //  TowerDigitizer->set_DoTowerDecal(true,"emcal_corr1_29.root",false);
  TowerDigitizer->set_DoTowerDecal(true,"emcal_newPatternCinco.root",false);

  // third parameter (doInverse) specifies if you want 
  //to instead apply the reciprocal  of the TowbyTow factors 
  // in the db file as a multiplicative factor 
  // here is an example of that
  //TowerDigitizer->set_DoTowerDecal(true,"emcal_newPatternCinco.root",true);

  // in actuality we do not actually suffer from digitization
  // effects on the entire pulse amplitude but rather sample
  //  ~12-14 times (pulse/pedestal itself about ~8 times)
  // and both the pedestal and pulse are extracted as continuous 
  // (or near continous) quantities.  In the near future the 
  // simple ("old fashioned") digitization scheme needs to implement
  // a full sim of the pulse extraction procedure. 
  // therefore for now when running in decal mode, as a temporary
  // more realistic approximation of this procedure, we change the 
  // energy response to continuous adc/energy values rather than digitized.
  // this behavior currently only occurs if running in the decal mode
  // ....
  // If you want to apply this non-digitizing part of the code
  // WITHOUT mod'ing the energy (ie w/o decal), call the DoDecal function without 
  // specifiying a filename as in the following example
  // 
  //  TowerDigitizer->set_DoTowerDecal(true, "",false);

  se->registerSubsystem(TowerDigitizer);
  

  RawTowerCalibration *TowerCalibration = new RawTowerCalibration("EmcRawTowerCalibration");
  TowerCalibration->Detector("CEMC");
  TowerCalibration->Verbosity(verbosity);
  



  if (G4CEMC::Cemc_spacal_configuration == PHG4CylinderGeom_Spacalv1::k1DProjectiveSpacal)
    {
    if (G4CEMC::TowerDigi == RawTowerDigitizer::kNo_digitization)
    {
      // just use sampling fraction set previously
      TowerCalibration->set_calib_const_GeV_ADC(1.0 / sampling_fraction);
    }
    else
    {
      TowerCalibration->set_calib_algorithm(RawTowerCalibration::kSimple_linear_calibration);
      TowerCalibration->set_calib_const_GeV_ADC(1. / photoelectron_per_GeV);
      TowerCalibration->set_pedstal_ADC(0);
    }
  }
  else if (G4CEMC::Cemc_spacal_configuration == PHG4CylinderGeom_Spacalv1::k2DProjectiveSpacal)
  {
    if (G4CEMC::TowerDigi == RawTowerDigitizer::kNo_digitization)
    {
      // just use sampling fraction set previously
      TowerCalibration->set_calib_const_GeV_ADC(1.0 / sampling_fraction);
    }
    else
    {
      
      if (!doSimple) 
	{


      //for tower by tower cal
      TowerCalibration->set_calib_algorithm(RawTowerCalibration::kTower_by_tower_calibration);
      TowerCalibration->GetCalibrationParameters().ReadFromFile("CEMC", "xml", 0, 0,
								string(getenv("CALIBRATIONROOT")) + string("/CEMC/TowerCalibCombinedParams_2020/"));  // calibration database
      TowerCalibration->set_variable_GeV_ADC(true);                                                                                                   //read GeV per ADC from calibrations file comment next line if true
      //    TowerCalibration->set_calib_const_GeV_ADC(1. / photoelectron_per_GeV / 0.9715);                                                             // overall energy scale based on 4-GeV photon simulations
      TowerCalibration->set_variable_pedestal(true);                                                                                                  //read pedestals from calibrations file comment next line if true
      //    TowerCalibration->set_pedstal_ADC(0);
      ///////////////////////////

	}
      else // dosimple
	{

	  TowerCalibration->set_calib_algorithm(RawTowerCalibration::kDbfile_tbt_gain_corr);
	  TowerCalibration->set_UseConditionsDB(false);

	  // Some standard db calo calibration files specification 
	  //  
	  // uses the db file accessor api and thus format 
	  // (format can change along with accessor internals, but user
	  // needs to know/give a readable format file).
	  // Decal tower by tow. factors in file are apply as a multiplicative 

	  //TowerCalibration->set_CalibrationFileName("emcal_corr1_29.root");
	  //TowerCalibration->set_CalibrationFileName("inv_emcal_corr_sec12bands.root");
	  TowerCalibration->set_CalibrationFileName("emcal_corr1_00.root");
	  //TowerCalibration->set_CalibrationFileName("emcal_newPatternCinco.root");

	  // since for this loop we avert the tower by tower improvements 
	  // in the kTower_by_tower xml-file based cemc calibration 
	  // which can be reimplemented in the new tower by tower dbfile format
	  // but for now we make a single overall reduction factor of 0.87
	  // that matches the same average calibration correction
	  // changing the following line to the one after it.  -JEF
	  //          TowerCalibration->set_calib_const_GeV_ADC(1. / photoelectron_per_GeV / 0.9715);                                                             // overall energy scale based on 4-GeV photon simulations
          TowerCalibration->set_calib_const_GeV_ADC(0.87 * 1./ photoelectron_per_GeV / 0.9715);                                                             // overall energy scale based on 4-GeV photon simulations
	  TowerCalibration->set_pedstal_ADC(0);
	  // note that in the TowerCalibration object, the pedestal subtraction is no
	  // longer applied, the above line simply follows suit with all the other 
	  // "calib_algorithms" on that object

	  ///////////////////////////
	}
      
    }
  }
  else
  {
    cout << "G4_CEmc_Spacal.C::CEMC_Towers - Fatal Error - unrecognized SPACAL configuration #"
         << G4CEMC::Cemc_spacal_configuration << ". Force exiting..." << endl;
    gSystem->Exit(-1);
    return;
  }
  se->registerSubsystem(TowerCalibration);

  return;
}

void CEMC_Clusters()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::CEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  if (G4CEMC::Cemc_clusterizer == G4CEMC::kCemcTemplateClusterizer)
  {
    RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->Verbosity(verbosity);
    ClusterBuilder->set_threshold_energy(0.030);  // This threshold should be the same as in CEMCprof_Thresh**.root file below
    std::string emc_prof = getenv("CALIBRATIONROOT");
    emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
    ClusterBuilder->LoadProfile(emc_prof);
    se->registerSubsystem(ClusterBuilder);
  }
  else if (G4CEMC::Cemc_clusterizer == G4CEMC::kCemcGraphClusterizer)
  {
    RawClusterBuilderGraph *ClusterBuilder = new RawClusterBuilderGraph("EmcRawClusterBuilderGraph");
    ClusterBuilder->Detector("CEMC");
    ClusterBuilder->Verbosity(verbosity);
    se->registerSubsystem(ClusterBuilder);
  }
  else
  {
    cout << "CEMC_Clusters - unknown clusterizer setting!" << endl;
    exit(1);
  }

  RawClusterPositionCorrection *clusterCorrection = new RawClusterPositionCorrection("CEMC");

  clusterCorrection->Get_eclus_CalibrationParameters().ReadFromFile("CEMC_RECALIB", "xml", 0, 0,
                                                                    //raw location
                                                                    string(getenv("CALIBRATIONROOT")) + string("/CEMC/PositionRecalibration_EMCal_9deg_tilt/"));

  clusterCorrection->Get_ecore_CalibrationParameters().ReadFromFile("CEMC_ECORE_RECALIB", "xml", 0, 0,
                                                                    //raw location
                                                                    string(getenv("CALIBRATIONROOT")) + string("/CEMC/PositionRecalibration_EMCal_9deg_tilt/"));

  clusterCorrection->Verbosity(verbosity);
  se->registerSubsystem(clusterCorrection);

  return;
}
void CEMC_Eval(const std::string &outputfile)
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::CEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();

  CaloEvaluator *eval = new CaloEvaluator("CEMCEVALUATOR", "CEMC", outputfile);
  eval->Verbosity(verbosity);
  se->registerSubsystem(eval);

  return;
}

void CEMC_QA()
{
  int verbosity = std::max(Enable::QA_VERBOSITY, Enable::CEMC_VERBOSITY);

  Fun4AllServer *se = Fun4AllServer::instance();
  QAG4SimulationCalorimeter *qa = new QAG4SimulationCalorimeter("CEMC");
  qa->Verbosity(verbosity);
  se->registerSubsystem(qa);

  return;
}

#endif
