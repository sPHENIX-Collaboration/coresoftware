// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EMCALSHOWERSHAPES_H
#define EMCALSHOWERSHAPES_H

#include "JetQADefs.h"

#include <fun4all/SubsysReco.h>

#include <limits>
#include <string>

class Fun4AllHistoManager;
class PHCompositeNode;
class RawCluster;
class RawTowerGeom;
class RawTowerGeomContainer;
class TH1;
class TH2;
class TowerInfoContainer;
class TriggerAnalyzer;

class EMCalShowerShapes : public SubsysReco
{
 public:
  EMCalShowerShapes(const std::string &modulename = "EMCalShowerShapes", const std::string &inputnode = "CLUSTERINFO_CEMC", const std::string &histtag = "");
  ~EMCalShowerShapes() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  //int ResetEvent(PHCompositeNode *topNode) override;
  //int EndRun(const int runnumber) override;
  //int End(PHCompositeNode *topNode) override;
  //int Reset(PHCompositeNode *topNode) override;
  void Print(const std::string &what = "ALL") const override;

  void SetTrgToSelect(const uint32_t trig = JetQADefs::GL1::MBDNSPhoton1)
  {
    m_doTrgSelect = true;
    m_trgToSelect = trig;
  }

  void SetHistTag(const std::string& tag)
  {
    m_histtag = tag;
  }

  void SetApplyMbdZvtxCut(const bool apply)
  {
    m_doMbdZvtxCut = apply;
  }

  void SetMbdZvtxMax(const float maxz)
  {
    m_mbdZvtxMax = maxz;
  }

  void SetApplyClusterEtaCut(const bool apply)
  {
    m_doClusterEtaCut = apply;
  }

  void SetApplyClusterETCut(const bool apply)
  {
    m_doClusterETCut = apply;
  }

  void SetClusterEtaMax(const float maxeta)
  {
    m_clusterEtaMax = maxeta;
  }

 private:
  struct ShowerShapeData
  {
    float e11 {std::numeric_limits<float>::quiet_NaN()};
    float e33 {std::numeric_limits<float>::quiet_NaN()};
    float e55 {std::numeric_limits<float>::quiet_NaN()};
    float e77 {std::numeric_limits<float>::quiet_NaN()};
    float e32 {std::numeric_limits<float>::quiet_NaN()};
    float e35 {std::numeric_limits<float>::quiet_NaN()};
    float weta {std::numeric_limits<float>::quiet_NaN()};
    float wphi {std::numeric_limits<float>::quiet_NaN()};
    float weta_cogx {std::numeric_limits<float>::quiet_NaN()};
    float wphi_cogx {std::numeric_limits<float>::quiet_NaN()};
    float detamax {std::numeric_limits<float>::quiet_NaN()};
    float dphimax {std::numeric_limits<float>::quiet_NaN()};
    float mean_time {std::numeric_limits<float>::quiet_NaN()};
    float iso04_emcal {std::numeric_limits<float>::quiet_NaN()};
  };

  bool LoadEMCalNodes(PHCompositeNode *topNode);
  float GetVertexZ(PHCompositeNode *topNode) const;
  bool CalculateShowerShapes(RawCluster* cluster, float cluster_eta, float cluster_phi, float cluster_et, float vertex_z, ShowerShapeData& data) const;
  double GetTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz) const;
  double DeltaR(double eta1, double phi1, double eta2, double phi2) const;
  float CalculateLayerET(float seed_eta, float seed_phi, float radius, TowerInfoContainer* towerContainer, RawTowerGeomContainer* geomContainer, float vertex_z) const;

  TriggerAnalyzer* m_analyzer {nullptr};
  Fun4AllHistoManager* m_manager {nullptr};
  TowerInfoContainer* m_emc_tower_container {nullptr};
  RawTowerGeomContainer* m_geomEM {nullptr};
  std::string m_modulename;
  std::string m_inputnode;
  std::string m_histtag;
  uint32_t m_trgToSelect;
  bool m_doTrgSelect;
  bool m_reportedMissingClusterNode {false};
  bool m_reportedMissingCaloNodes {false};
  float m_shape_min_tower_E {0.070F};
  bool m_doMbdZvtxCut {true};
  float m_mbdZvtxMax {60.0F};
  bool m_doClusterEtaCut {true};
  float m_clusterEtaMax {0.7F};
  bool m_doClusterETCut {true};
  float m_clusterETMin {5.0F};

  TH1* h_cluster_et {nullptr};
  TH1* h_e11oe33 {nullptr};
  TH1* h_e33oe55 {nullptr};
  TH1* h_e55oe77 {nullptr};
  TH1* h_e32oe35 {nullptr};
  TH1* h_weta {nullptr};
  TH1* h_wphi {nullptr};
  TH1* h_weta_cogx {nullptr};
  TH1* h_wphi_cogx {nullptr};
  TH1* h_detamax {nullptr};
  TH1* h_dphimax {nullptr};
  TH1* h_mean_time {nullptr};
  TH1* h_iso04_emcal {nullptr};
  TH2* h_weta_vs_et {nullptr};
  TH2* h_wphi_vs_et {nullptr};
};

#endif
