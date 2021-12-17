#ifndef CALORECO_RAWTOWERCALIBRATION_H
#define CALORECO_RAWTOWERCALIBRATION_H

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameters.h>

#include <array>
#include <cmath>
#include <iostream>
#include <string>

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
  ~RawTowerCalibration() override
  {
  }

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void Detector(const std::string &d)
  {
    m_Detector = d;
    m_TowerCalibParams.set_name(d);
  }

  void CalibFile(const std::string &f)
  {
    m_CalibrationFileName = f;
  }

  void TowerType(const int type)
  {
    m_TowerType = type;
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

  enu_calib_algorithm get_calib_algorithm() const
  {
    return m_CalibAlgorithmEnum;
  }

  void set_calib_algorithm(enu_calib_algorithm calibAlgorithm)
  {
    m_CalibAlgorithmEnum = calibAlgorithm;
  }

  double get_calib_const_GeV_ADC() const
  {
    return m_CalibConst_GeV_per_ADC;
  }

  void set_calib_const_GeV_ADC(double calibConstGeVAdc)
  {
    m_CalibConst_GeV_per_ADC = calibConstGeVAdc;
  }

  void set_variable_GeV_ADC(const bool value)
  {
    m_GeV_per_ADC_FromFileFlag = value;
  }

  std::string get_calib_tower_node_prefix() const
  {
    return m_CalibTowerNodePrefix;
  }

  void set_calib_tower_node_prefix(const std::string &calibTowerNodePrefix)
  {
    m_CalibTowerNodePrefix = calibTowerNodePrefix;
  }

  double get_pedstal_ADC() const
  {
    return m_PedestalADC;
  }

  void set_pedstal_ADC(const double pedstalAdc)
  {
    m_PedestalADC = pedstalAdc;
  }

  void set_variable_pedestal(const bool value)
  {
    m_PedestalFromFileFlag = value;
  }

  std::string get_raw_tower_node_prefix() const
  {
    return m_RawTowerNodePrefix;
  }

  void set_raw_tower_node_prefix(const std::string &rawTowerNodePrefix)
  {
    m_RawTowerNodePrefix = rawTowerNodePrefix;
  }

  //! Get the parameters for update. Useful fields are listed in SetDefaultParameters();
  PHParameters &GetCalibrationParameters()
  {
    return m_TowerCalibParams;
  }

 protected:
  void CreateNodes(PHCompositeNode *topNode);

  RawTowerContainer *m_CalibTowerContainer = nullptr;
  RawTowerContainer *m_RawTowerContainer = nullptr;
  RawTowerGeomContainer *m_RawTowerGeomContainer = nullptr;

  enu_calib_algorithm m_CalibAlgorithmEnum = kNo_calibration;

  //! tower type to act on
  int m_TowerType = -1;

  //! pedstal in unit of ADC
  double m_PedestalADC = NAN;

  //! calibration constant in unit of GeV per ADC
  double m_CalibConst_GeV_per_ADC = NAN;

//  double m_RecalArray[24][64] = {{0.0}};
  std::array< std::array<double,64>, 24> m_RecalArray;
  //! pedestal from file
  bool m_PedestalFromFileFlag = false;

  //! GeV per ADC from file
  bool m_GeV_per_ADC_FromFileFlag = false;


  std::string m_Detector;
  std::string m_CalibrationFileName;

  std::string m_CalibTowerNodePrefix = "CALIB";
  std::string m_RawTowerNodePrefix = "RAW";

  //! Tower by tower calibration parameters
  PHParameters m_TowerCalibParams;

};

#endif
