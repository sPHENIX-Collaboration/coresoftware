#ifndef QA_QAG4SIMULATIONMVTX_H
#define QA_QAG4SIMULATIONMVTX_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <memory>
#include <string>
#include <set>
#include <utility>

class PHG4Hit;
class PHG4HitContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;

/// \class QAG4SimulationMvtx
class QAG4SimulationMvtx : public SubsysReco
{

  public:

  /// constructor
  QAG4SimulationMvtx(const std::string & name = "QAG4SimulationMvtx");

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  private:

  /// common prefix for QA histograms
  std::string get_histo_prefix() const;

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// evaluate clusters
  void evaluate_clusters();

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;

  /// cluster map
  TrkrClusterContainer* m_cluster_map = nullptr;

  /// clusters to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

  /// hit to g4hit association
  TrkrHitTruthAssoc* m_hit_truth_map = nullptr;

  /// g4 hits
  PHG4HitContainer* m_g4hits_mvtx = nullptr;

  /// list of relevant layers
  /* it is filled at Init stage. It should not change for the full run */
  std::set<int> m_layers;

};

#endif
