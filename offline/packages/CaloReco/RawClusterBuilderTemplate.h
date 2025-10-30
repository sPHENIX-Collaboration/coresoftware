#ifndef CALORECO_RAWCLUSTERBUILDERTEMPLATE_H
#define CALORECO_RAWCLUSTERBUILDERTEMPLATE_H

#include <fun4all/SubsysReco.h>

#include <globalvertex/GlobalVertex.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;
class RawClusterv2;
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
  void PrintCylGeom(RawTowerGeomContainer* towergeom, const std::string& fname) const;
  void SetProfileProb(bool pprob) { bProfProb = pprob; }
  void SetProbNoiseParam(float rn) { fProbNoiseParam = rn; }
    
  void set_threshold_energy(const float e) { _min_tower_e = e; }
  void set_peakthreshold_energy(const float e) { _min_peak_e = e; }
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
  {  // 0 use global vtx, 1 only bbcout bbczvtx, 2 use NO zvtx[set to 0], 3 use MC truth vertex
    m_UseAltZVertex = useAltZMode;
  }

  void set_UseCorrectPosition(const bool useCorrectPosition);
  
  void set_UseCorrectShowerDepth(const bool useCorrectShowerDepth);

  void set_UseDetailedGeometry(const bool useDetailedGeometry);
    
  void WriteClusterV2(bool enable) { m_writeClusterV2 = enable; }

  void setOutputClusterNodeName(const std::string& inpNodenm)
  {
    m_outputnodename = inpNodenm;
  }

  void setTowerGeomNodeName(const std::string& name)
  {
    m_TowerGeomNodeName = name;
  }

  // !!! note :  next fn NOT implemented for RawTowers
  // only TowerInfo  mode
  void setInputTowerNodeName(const std::string& inpNodenm)
  {
    m_inputnodename = inpNodenm;
  }

  void set_min_cluster_E_saved(float min_cluster_E) { m_min_cluster_e = min_cluster_E; }

  void set_vertex_type(GlobalVertex::VTXTYPE type)
  {
    m_vertex_type = type;
  }

  void setSubclusterSplitting(bool doSubClusterSplitting)
  {
    m_subclustersplitting = doSubClusterSplitting;
  }

  

 private:
  void CreateNodes(PHCompositeNode* topNode);
  static bool Cell2Abs(RawTowerGeomContainer* towergeom, float phiC, float etaC, float& phi, float& eta);
  bool IsAcceptableTower(TowerInfo *tower) const;
  bool IsAcceptableTower(RawTower *tower) const;

  RawClusterContainer* _clusters{nullptr};
  //  BEmcProfile *_emcprof;

  BEmcRec* bemc{nullptr};
  float fEnergyNorm{1.};

  float _min_tower_e{0.020};
  float _min_peak_e{0.200};
  int chkenergyconservation{0};

  std::string detector;
  std::string ClusterNodeName;

  int BINX0{0};
  int NBINX{0};
  int BINY0{0};
  int NBINY{0};

  bool bPrintGeom{false};
  bool bProfProb{false};
  float fProbNoiseParam{0.04};

  GlobalVertex::VTXTYPE m_vertex_type{GlobalVertex::MBD};

  int m_UseTowerInfo{0};  // 0 only old tower, 1 only new (TowerInfo based),

  bool m_do_tower_selection{true};

  std::string m_towerInfo_nodename;


  bool m_UseDetailedGeometry {true};
  // Use a more detailed calorimeter geometry (default)
  // Only available for CEMC

  int m_UseAltZVertex{1};
  // 0 - use GlobalVtxMap
  // 1 - use BbcReco ZVtx
  // 2 - use NO zvertex (zvtx = 0)
  // 3 - use truth MC zvertex

  float m_min_cluster_e{0.0};

  bool m_subclustersplitting{true};
  bool m_writeClusterV2{false};  // default off: use RawClusterv1

  std::string m_inputnodename;
  std::string m_outputnodename;
  std::string m_TowerGeomNodeName;
};

#endif /* RawClusterBuilderTemplate_H__ */
