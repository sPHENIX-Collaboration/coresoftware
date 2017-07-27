#ifndef __RAWCLUSTERPOSITIONCORRECTION_H__
#define __RAWCLUSTERPOSITIONCORRECTION_H__

#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4Parameters.h>
#include <phool/PHObject.h>
#include <string>

class RawClusterContainer;
class RawCluster;
class RawTowerContainer;
class RawTower;
class RawClusterv1;

class RawClusterPositionCorrection : public SubsysReco
{
 public:
  explicit RawClusterPositionCorrection(const std::string &name);

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void CreateNodeTree(PHCompositeNode *topNode);

  const PHG4Parameters &GetCalibrationParameters() const
  {
    return _calib_params;
  }
  PHG4Parameters &GetCalibrationParameters()
  {
    return _calib_params;
  }

  void SetCalibrationParameters(const PHG4Parameters &calib_params)
  {
    _calib_params = calib_params;
  }

 private:
  PHG4Parameters _calib_params;

  void SetDefaultParameters(PHG4Parameters &param);
  RawClusterContainer *_recalib_clusters;

  std::string _det_name;

  int bins;
  std::vector<float> binvals;
  std::vector<std::vector<double> > calib_constants;
};

#endif  // __RAWCLUSTERPOSITIONCORRECTION_H__
