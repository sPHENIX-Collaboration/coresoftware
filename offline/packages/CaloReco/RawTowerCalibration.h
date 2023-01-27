#ifndef CALORECO_RAWTOWERCALIBRATION_H
#define CALORECO_RAWTOWERCALIBRATION_H

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameters.h>

#include <cmath>
#include <iostream>
#include <string>

class CaloCalibSimpleCorrFile;
class PHCompositeNode;
class RawTowerContainer;
class TowerInfoContainerv1;
class RawTowerGeomContainer;

//! calibrate ADC value to measured energy deposition in calorimeter towers
//! default input DST node is TOWER_RAW_DETECTOR
//! default output DST node is TOWER_CALIB_DETECTOR
class RawTowerCalibration : public SubsysReco
{
 public:
  RawTowerCalibration(const std::string &name = "RawTowerCalibration");
  ~RawTowerCalibration() override
  {
  }

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void
  Detector(const std::string &d)
  {
    detector = d;
    _tower_calib_params.set_name(d);
  }
  void
  TowerType(const int type)
  {
    _tower_type = type;
  }

  enum enu_calib_algorithm
  {
    //! directly pass the energy of raw tower to calibrated tower. Zero suppression is applied
    kNo_calibration = 0,

    //! simple calibration with pedstal subtraction and a global energy scale (sampling fraction) correction
    kSimple_linear_calibration = 1,

    //! input calibration file for tower by tower calibration. Use GetCalibrationParameters() to set the calibration parameters
    kTower_by_tower_calibration = 2,

    // use conditions DB file/wrapper (non-xml) file for most gain tracing correction factors
    kDbfile_tbt_gain_corr = 3
  };
  enum ProcessTowerType
  {
    kRawTowerOnly= 0,
    kTowerInfoOnly = 1,
    kBothTowers =2
  };

  enu_calib_algorithm
  get_calib_algorithm() const
  {
    return _calib_algorithm;
  }

  void
  set_calib_algorithm(enu_calib_algorithm calibAlgorithm)
  {
    _calib_algorithm = calibAlgorithm;
  }

  double
  get_calib_const_GeV_ADC() const
  {
    return _calib_const_GeV_ADC;
  }

  void
  set_calib_const_GeV_ADC(double calibConstGeVAdc)
  {
    _calib_const_GeV_ADC = calibConstGeVAdc;
  }

  void
  set_variable_GeV_ADC(const bool value)
  {
    _GeV_ADC_file = value;
  }

  std::string
  get_calib_tower_node_prefix() const
  {
    return _calib_tower_node_prefix;
  }

  void
  set_calib_tower_node_prefix(const std::string &calibTowerNodePrefix)
  {
    _calib_tower_node_prefix = calibTowerNodePrefix;
  }

  double
  get_pedstal_ADC() const
  {
    return _pedstal_ADC;
  }

  void
  set_pedstal_ADC(double pedstalAdc)
  {
    _pedstal_ADC = pedstalAdc;
  }

  void
  set_variable_pedestal(const bool value)
  {
    _pedestal_file = value;
  }

  std::string
  get_raw_tower_node_prefix() const
  {
    return _raw_tower_node_prefix;
  }

  void
  set_raw_tower_node_prefix(const std::string &rawTowerNodePrefix)
  {
    _raw_tower_node_prefix = rawTowerNodePrefix;
  }

  void
  set_zero_suppression_GeV(double)
  {
    std::cout << "RawTowerCalibration::set_zero_suppression_GeV is deprecated!" << std::endl
              << "  See discussion at https://github.com/sPHENIX-Collaboration/coresoftware/pull/867" << std::endl
              << std::endl;
  }

  //! Get the parameters for update. Useful fields are listed in SetDefaultParameters();
  PHParameters &
  GetCalibrationParameters()
  {
    return _tower_calib_params;
  }

  void set_CalibrationFileName(const char *inCalFname)
  {
    m_CalibrationFileName = inCalFname;
  }
  void set_UseConditionsDB(const bool setUseCondDB)
  {
    m_UseConditionsDB = setUseCondDB;
  }


  void set_towerinfo(RawTowerCalibration::ProcessTowerType UseTowerInfo )
  {
    m_UseTowerInfo = UseTowerInfo;
  }


 protected:
  void
  CreateNodes(PHCompositeNode *topNode);

  enu_calib_algorithm _calib_algorithm;

  RawTowerContainer *_calib_towers = nullptr;
  RawTowerContainer *_raw_towers = nullptr;

  TowerInfoContainerv1 *_calib_towerinfos = nullptr;
  TowerInfoContainerv1 *_raw_towerinfos = nullptr;


  RawTowerGeomContainer *rawtowergeom = nullptr;

  std::string detector;
  std::string RawTowerNodeName;
  std::string RawTowerInfoNodeName;
  std::string CaliTowerNodeName;
  std::string CaliTowerInfoNodeName;
  std::string TowerGeomNodeName;

  std::string _calib_tower_node_prefix;
  std::string _raw_tower_node_prefix;

  //! pedstal in unit of ADC
  double _pedstal_ADC = NAN;

  //! pedestal from file
  bool _pedestal_file = false;

  //! calibration constant in unit of GeV per ADC
  double _calib_const_GeV_ADC = NAN;

  //! GeV per ADC from file
  bool _GeV_ADC_file = false;

  //! tower type to act on
  int _tower_type = -1;

  //! Tower by tower calibration parameters
  PHParameters _tower_calib_params;

  std::string m_CalibrationFileName;
  bool m_UseConditionsDB = false;

  CaloCalibSimpleCorrFile *_cal_dbfile = nullptr;
  RawTowerCalibration::ProcessTowerType m_UseTowerInfo = RawTowerCalibration::ProcessTowerType::kBothTowers;  // 0 just produce RawTowers, 1 just produce TowerInfo objects, and 2 produce both


};

#endif
