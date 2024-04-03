#ifndef CALORECO_RAWCLUSTERPOSITIONCORRECTION_H
#define CALORECO_RAWCLUSTERPOSITIONCORRECTION_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class PHCompositeNode;
class RawClusterContainer;
class CDBHistos;
class CDBInterface;
class CDBTTree;
class TH1;
class TH2;

class RawClusterPositionCorrection : public SubsysReco
{
 public:
  explicit RawClusterPositionCorrection(const std::string &name);
  ~RawClusterPositionCorrection() override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void CreateNodeTree(PHCompositeNode *topNode);

  // const PHParameters &Get_eclus_CalibrationParameters() const
  // {
  //   return _eclus_calib_params;
  // }
  // PHParameters &Get_eclus_CalibrationParameters()
  // {
  //   return _eclus_calib_params;
  // }
  // const PHParameters &Get_eore_CalibrationParameters() const
  // {
  //   return _ecore_calib_params;
  // }
  // PHParameters &Get_ecore_CalibrationParameters()
  // {
  //   return _ecore_calib_params;
  // }

  // void Set_eclus_CalibrationParameters(const PHParameters &calib_params)
  // {
  //   _eclus_calib_params = calib_params;
  // }
  // void Set_ecore_CalibrationParameters(const PHParameters &calib_params)
  // {
  //   _ecore_calib_params = calib_params;
  // }

  void set_UseTowerInfo(const int useMode)
  {  // 0 only old tower, 1 only new (TowerInfo based),
    m_UseTowerInfo = useMode;
  }

 private:
  // PHParameters _eclus_calib_params;
  // PHParameters _ecore_calib_params;
  // void SetDefaultParameters(PHParameters &param);
  RawClusterContainer *_recalib_clusters{};

  std::string _det_name;

  // std::vector<float> binvals;
  // std::vector<std::vector<double> > eclus_calib_constants;
  // std::vector<std::vector<double> > ecore_calib_constants;

  // key: phibin, etabin
  std::vector<std::vector<float>> calib_constants_north;
  std::vector<std::vector<float>> calib_constants_north_ecore;
  std::vector<std::vector<float>> calib_constants_south;
  std::vector<std::vector<float>> calib_constants_south_ecore;

  int m_UseTowerInfo {0};  // 0 only old tower, 1 only new (TowerInfo based),

  int bins_eta;
  int bins_phi;
  int iEvent;

  TH2* h2NorthSector{nullptr};
  TH2* h2SouthSector{nullptr};
  TH1* pdcCorrFlat{nullptr};
  
  CDBTTree *cdbttree{nullptr};
  CDBHistos *cdbHisto{nullptr};
};

#endif
