#ifndef __CaloCalibrationF__
#define __CaloCalibrationF__

//* Unpacks raw HCAL PRDF files *//
//Abhisek Sen

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <string>
#include <phparameter/PHParameters.h>

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

  //! Get the parameters for readonly
  const PHParameters &
  GetCalibrationParameters() const
  {
    return _calib_params;
  }

  //! Get the parameters for update. Useful fields are listed in SetDefaultParameters();
  PHParameters &
  GetCalibrationParameters()
  {
    return _calib_params;
  }

  //! Overwrite the parameter. Useful fields are listed in SetDefaultParameters();
  void
  SetCalibrationParameters(const PHParameters & calib_params)
  {
    _calib_params = calib_params;
  }

private:

  RawTowerContainer* _calib_towers;
  RawTowerContainer* _raw_towers;

  std::string detector;
  std::string RawTowerNodeName;
  std::string CaliTowerNodeName;

  std::string _calib_tower_node_prefix;
  std::string _raw_tower_node_prefix;

  PHParameters _calib_params;

  //! load the default parameter to param
  void
  SetDefaultParameters(PHParameters & param);

};

#endif //**CaloCalibrationF**//
