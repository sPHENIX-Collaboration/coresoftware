#ifndef CALOVTXALGOCaloZ_H
#define CALOVTXALGOCaloZ_H

#include "CaloVtxAlgo.h"
#include <string>

class PHCompositeNode;

// Wraps the trained CaloZ (VertexCaloZ.h / vertex_mlp_weights.root) as a
// CaloVtxAlgo so it can be registered alongside other algorithms in
// CaloVtxReco and compared directly.
class CaloVtxAlgoCaloZ : public CaloVtxAlgo
{
 public:
  explicit CaloVtxAlgoCaloZ();
  ~CaloVtxAlgoCaloZ() override = default;

  int Init(PHCompositeNode *topNode) override;
  int CalculateVertex(PHCompositeNode *topNode, float &zvtx) override;
  std::string Name() const override { return "CaloZ"; }
  
  float get_energy_cut() { return m_energy_cut; }
  void set_energy_cut(float new_energy) { m_energy_cut = new_energy; }

 private:

  float m_energy_cut{0.1};
};

#endif
