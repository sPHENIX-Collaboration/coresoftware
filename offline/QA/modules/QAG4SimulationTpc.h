#ifndef QA_QAG4SIMULATIONTPC_H
#define QA_QAG4SIMULATIONTPC_H

#include <g4eval/SvtxEvalStack.h>  // for SvtxEvalStack

#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <map>
#include <memory>
#include <set>
#include <string>

class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHG4TruthInfoContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class SvtxEvalStack;
class ActsGeometry;

/// \class QAG4SimulationTpc
class QAG4SimulationTpc : public SubsysReco
{
 public:
  /// constructor
  QAG4SimulationTpc(const std::string& name = "QAG4SimulationTpc");

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  void set_cluster_version(int value) { m_cluster_version = value; }

 private:
  /// common prefix for QA histograms
  std::string get_histo_prefix() const;

  std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;
  PHG4TruthInfoContainer* m_truthContainer;

  /// load nodes
  int load_nodes(PHCompositeNode*);

  /// evaluate clusters
  void evaluate_clusters();

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits(TrkrDefs::cluskey) const;

  /// true if histograms are initialized
  bool m_initialized = false;

  /// Acts tracking geometry for surface lookup
  ActsGeometry* m_tGeometry = nullptr;

  /// cluster map
  TrkrClusterContainer* m_cluster_map = nullptr;

  /// clusters to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

  /// hit to g4hit association
  TrkrHitTruthAssoc* m_hit_truth_map = nullptr;

  /// g4 hits
  PHG4HitContainer* m_g4hits_tpc = nullptr;

  /// list of relevant layers
  /* it is filled at Init stage. It should not change for the full run */
  std::set<int> m_layers;
  std::multimap<int, int> m_layer_region_map;
  ClusterErrorPara _ClusErrPara;
  int m_cluster_version = 4;
};

#endif
