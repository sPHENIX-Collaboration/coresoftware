#ifndef CALOVTXALGOMLP_H
#define CALOVTXALGOMLP_H

#include "CaloVtxAlgo.h"

#include <string>

class PHCompositeNode;

// Wraps the trained MLP (VertexMLP.h / vertex_mlp_weights.root) as a
// CaloVtxAlgo so it can be registered alongside other algorithms in
// CaloVtxReco and compared directly.
class CaloVtxAlgoMLP : public CaloVtxAlgo
{
 public:
  explicit CaloVtxAlgoMLP() = default;
  ~CaloVtxAlgoMLP() override = default;

  int Init(PHCompositeNode *topNode) override;
  int CalculateVertex(PHCompositeNode *topNode, float &zvtx) override;
  std::string Name() const override { return "MLP"; }
  VertexDefs::CALOALGO Algo() const override { return VertexDefs::CALOALGO::JETMLP; }
  
  void setWeightsFile(std::string file) { m_weightsFile = file; }
  void setJetNode(std::string node) { m_jet_node = node; }
private:

  std::string m_jet_node{"Antikt_TowerInfo_vtx_none_r06"};
  std::string m_weightsFile{"vertex_mlp_weights.root"};

  float m_radius_EM{0};
  float m_radius_IH{0};
  float m_radius_OH{0};
};

#endif
