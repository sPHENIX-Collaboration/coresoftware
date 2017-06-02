#ifndef __RAWCLUSTERRECALIBRATOR_H__
#define __RAWCLUSTERRECALIBRATOR_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <string>
#include <g4detectors/PHG4Parameters.h>

class RawClusterContainer;
class RawCluster;
class RawTowerContainer;
class RawTower;
class RawClusterv1;


class RawClusterRecalibrator: public SubsysReco
{
public:
  RawClusterRecalibrator(const std::string& name);
  
  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void CreateNodeTree(PHCompositeNode *topNode);

  
  const PHG4Parameters & GetCalibrationParameters() const
  {
    return _calib_params;
  }
  PHG4Parameters & GetCalibrationParameters()
  {
    return _calib_params;
  }
  
  void SetCalibrationParameters(const PHG4Parameters & calib_params)
  {
    _calib_params = calib_params;
  }
  


private:

  PHG4Parameters _calib_params;

  void SetDefaultParameters(PHG4Parameters & param);
  RawClusterContainer *_recalib_clusters;

  int bins;
  std::vector<float> binvals;

};

#endif // __RAWCLUSTERRECLIBRATOR_H__
