#ifndef CALOVTXALGOJetSkew_H
#define CALOVTXALGOJetSkew_H

#include "CaloVtxAlgo.h"
#include <calobase/RawTowerDefs.h>
#include <string>
#include "TMath.h"
#include <cmath>
#include <memory>


class PHCompositeNode;
class RawTowerGeomContainer;
class TowerInfoContainer;

// Wraps the trained JetSkew (VertexJetSkew.h / vertex_mlp_weights.root) as a
// CaloVtxAlgo so it can be registered alongside other algorithms in
// CaloVtxReco and compared directly.
class CaloVtxAlgoJetSkew : public CaloVtxAlgo
{
 public:
  explicit CaloVtxAlgoJetSkew() = default;
  ~CaloVtxAlgoJetSkew() override = default;

  int Init(PHCompositeNode *topNode) override;
  int CalculateVertex(PHCompositeNode *topNode, float &zvtx) override;
  std::string Name() const override { return "JetSkew"; }
  VertexDefs::CALOALGO Algo() const override { return VertexDefs::CALOALGO::JETSKEW; }
  
  void setJetNode(std::string jetnode) { m_jetnodename = jetnode; }
  float get_jet_threshold() { return m_jet_threshold; }
  void set_jet_threshold(float new_thresh) { m_jet_threshold = new_thresh; }
  float get_calib_factor() { return m_calib_factor; }
  void set_calib_factor(float new_calib) { m_calib_factor = new_calib; }
  float get_energy_cut() { return m_energy_cut; }
  void set_energy_cut(float new_energy) { m_energy_cut = new_energy; }

 private:

  float new_eta(int channel, TowerInfoContainer *towerset, RawTowerGeomContainer *geom, RawTowerDefs::CalorimeterId caloID, float testz);
  
  std::string m_jetnodename;
  float m_jet_threshold{15};
  float m_calib_factor{1.406};
  float m_energy_cut{0.1};
  float m_radius_EM{std::numeric_limits<float>::quiet_NaN()};
  float m_radius_OH{std::numeric_limits<float>::quiet_NaN()};

  TowerInfoContainer *towers[3] = {0};
  RawTowerGeomContainer *geom[3] = {0};


};

#endif
