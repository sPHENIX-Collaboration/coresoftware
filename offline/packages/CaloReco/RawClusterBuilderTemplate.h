#ifndef CALORECO_RAWCLUSTERBUILDERTEMPLATE_H
#define CALORECO_RAWCLUSTERBUILDERTEMPLATE_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;
class RawTowerGeomContainer;
class BEmcRec;
class TowerInfo;
class RawTower;

class RawClusterBuilderTemplate : public SubsysReco
{
 public:
  explicit RawClusterBuilderTemplate(const std::string& name = "RawClusterBuilderTemplate");
  ~RawClusterBuilderTemplate() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  void Detector(const std::string& d);

  void SetCylindricalGeometry();
  void SetPlanarGeometry();
  void PrintGeometry() { bPrintGeom = true; }  // Prints it at InitRun time
  void PrintCylGeom(RawTowerGeomContainer* towergeom, const std::string& fname);
  void SetProfileProb(bool pprob) { bProfProb = pprob; }
  void SetProbNoiseParam(float rn) { fProbNoiseParam = rn; }

  void set_threshold_energy(const float e) { _min_tower_e = e; }
  void setEnergyNorm(const float norm) { fEnergyNorm = norm; }
  void checkenergy(const int i = 1) { chkenergyconservation = i; }
  void LoadProfile(const std::string& fname);

  void set_UseTowerInfo(const int useMode)
  {  // 0 only old tower, 1 only new (TowerInfo based),
    m_UseTowerInfo = useMode;
  }

  void set_ApplyTowerSelection(bool b)
  {
    m_do_tower_selection = b;
  }

  void set_UseAltZVertex(const int useAltZMode)
  {  // 0 use global vtx, 1 only bbcout bbczvtx , 2 use NO zvtx[set to 0]
    m_UseAltZVertex = useAltZMode;
  }

  void setOutputClusterNodeName(const std::string& inpNodenm)
  {
    m_outputnodename = inpNodenm;
  }

  // !!! note :  next fn NOT implemented for RawTowers
  // only TowerInfo  mode
  void setInputTowerNodeName(const std::string& inpNodenm)
  {
    m_inputnodename = inpNodenm;
  }

 private:
  void CreateNodes(PHCompositeNode* topNode);
  bool Cell2Abs(RawTowerGeomContainer* towergeom, float phiC, float etaC, float& phi, float& eta);
  bool IsAcceptableTower(TowerInfo *tower);
  bool IsAcceptableTower(RawTower *tower);

  RawClusterContainer* _clusters = nullptr;
  //  BEmcProfile *_emcprof;

  BEmcRec* bemc = nullptr;
  float fEnergyNorm = 1.;

  float _min_tower_e = 0.020;
  int chkenergyconservation = 0;

  std::string detector;
  std::string ClusterNodeName;

  int BINX0 = 0;
  int NBINX = 0;
  int BINY0 = 0;
  int NBINY = 0;

  bool bPrintGeom = false;
  bool bProfProb = false;
  float fProbNoiseParam = 0.04;

  int m_UseTowerInfo = 0;  // 0 only old tower, 1 only new (TowerInfo based),

  bool m_do_tower_selection = true;

  std::string m_towerInfo_nodename;

  int m_UseAltZVertex = 2;
  // 0 - use GlobalVtxMap
  // 1 - use BbcReco ZVtx
  // 2 - use NO zvertex (zvtx = 0)

  std::string m_inputnodename;
  std::string m_outputnodename;
};

#endif /* RawClusterBuilderTemplate_H__ */
