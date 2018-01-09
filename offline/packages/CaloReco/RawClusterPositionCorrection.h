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

  const PHG4Parameters &Get_eclus_CalibrationParameters() const
  {
    return _eclus_calib_params;
  }
  PHG4Parameters &Get_eclus_CalibrationParameters()
  {
    return _eclus_calib_params;
  }
  const PHG4Parameters &Get_eore_CalibrationParameters() const
  {
    return _ecore_calib_params;
  }
  PHG4Parameters &Get_ecore_CalibrationParameters()
  {
    return _ecore_calib_params;
  }

  void Set_eclus_CalibrationParameters(const PHG4Parameters &calib_params)
  {
    _eclus_calib_params = calib_params;
  }
  void Set_ecore_CalibrationParameters(const PHG4Parameters &calib_params)
  {
    _ecore_calib_params = calib_params;
  }

 private:
  PHG4Parameters _eclus_calib_params;
  PHG4Parameters _ecore_calib_params;
  void SetDefaultParameters(PHG4Parameters &param);
  RawClusterContainer *_recalib_clusters;

  std::string _det_name;

  int bins;
  std::vector<float> binvals;
  std::vector<std::vector<double> > eclus_calib_constants;
  std::vector<std::vector<double> > ecore_calib_constants;
};

#endif  // __RAWCLUSTERPOSITIONCORRECTION_H__
