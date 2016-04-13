#ifndef __CaloCalibrationF__
#define __CaloCalibrationF__

//* Unpacks raw HCAL PRDF files *//
//Abhisek Sen

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <string>

class RawTowerContainer;

class CaloCalibration : public SubsysReco
{
public:
  CaloCalibration(const std::string& name);

  int
  Init(PHCompositeNode *topNode);

  int
  InitRun(PHCompositeNode *topNode);

  int
  process_event(PHCompositeNode *topNode);

  int
  End(PHCompositeNode *topNode);

  void
  CreateNodeTree(PHCompositeNode *topNode);

  std::string
  get_calib_tower_node_prefix() const
  {
    return _calib_tower_node_prefix;
  }

  void
  set_calib_tower_node_prefix(std::string calibTowerNodePrefix)
  {
    _calib_tower_node_prefix = calibTowerNodePrefix;
  }

  std::string
  get_raw_tower_node_prefix() const
  {
    return _raw_tower_node_prefix;
  }

  void
  set_raw_tower_node_prefix(std::string rawTowerNodePrefix)
  {
    _raw_tower_node_prefix = rawTowerNodePrefix;
  }

  double
  get_calib_const_scale() const
  {
    return _calib_const_scale;
  }

  void
  set_pedstal_ADC(double calib_const_scale)
  {
    _calib_const_scale = calib_const_scale;
  }
private:

  //! additional scale for the calibration constant
  double _calib_const_scale;

  RawTowerContainer* _calib_towers;
  RawTowerContainer* _raw_towers;

  std::string detector;
  std::string RawTowerNodeName;
  std::string CaliTowerNodeName;

  std::string _calib_tower_node_prefix;
  std::string _raw_tower_node_prefix;


};

#endif //**CaloCalibrationF**//
