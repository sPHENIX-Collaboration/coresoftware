#ifndef CALORECO_RAWTOWERCALIBRATION_H
#define CALORECO_RAWTOWERCALIBRATION_H

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameters.h>

#include <string>
#include <iostream>

class PHCompositeNode;
class RawTowerContainer;
class RawTowerGeomContainer;

//! calibrate ADC value to measured energy deposition in calorimeter towers
//! default input DST node is TOWER_RAW_DETECTOR
//! default output DST node is TOWER_CALIB_DETECTOR
class RawTowerCalibration : public SubsysReco
{
 public:
  RawTowerCalibration(const std::string &name = "RawTowerCalibration");
  virtual ~RawTowerCalibration()
  {
  }

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
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
    kTower_by_tower_calibration = 2
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
    std::cout <<"RawTowerCalibration::set_zero_suppression_GeV is deprecated!"<<std::endl
        <<"  See discussion at https://github.com/sPHENIX-Collaboration/coresoftware/pull/867"<<std::endl<<std::endl;
  }

  //! Get the parameters for update. Useful fields are listed in SetDefaultParameters();
  PHParameters &
  GetCalibrationParameters()
  {
    return _tower_calib_params;
  }

 protected:
  void
  CreateNodes(PHCompositeNode *topNode);

  enu_calib_algorithm _calib_algorithm;

  RawTowerContainer *_calib_towers;
  RawTowerContainer *_raw_towers;
  RawTowerGeomContainer *rawtowergeom;

  std::string detector;
  std::string RawTowerNodeName;
  std::string CaliTowerNodeName;
  std::string TowerGeomNodeName;

  std::string _calib_tower_node_prefix;
  std::string _raw_tower_node_prefix;

  //! pedstal in unit of ADC
  double _pedstal_ADC;

  //! calibration constant in unit of GeV per ADC
  double _calib_const_GeV_ADC;

  //! tower type to act on
  int _tower_type;

  //! Tower by tower calibration parameters
  PHParameters _tower_calib_params;
};

#endif
