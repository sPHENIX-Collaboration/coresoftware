#ifndef QA_QAG4SIMULATIONMICROMEGAS_H
#define QA_QAG4SIMULATIONMICROMEGAS_H

#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class ActsGeometry;
class PHCompositeNode;
class PHG4CylinderGeomContainer;
class PHG4Hit;
class PHG4Particle;
class PHG4HitContainer;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class PHG4TruthInfoContainer;
class SvtxEvalStack;

/// \class QAG4SimulationMicromegas
class QAG4SimulationMicromegas : public SubsysReco
{
 public:
  /// constructor
  QAG4SimulationMicromegas(const std::string& name = "QAG4SimulationMicromegas");

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  void set_cluster_version(int value) { m_cluster_version = value; }

 private:
  /// common prefix for QA histograms
  std::string get_histo_prefix() const;
  std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;
  /// load nodes
  int load_nodes(PHCompositeNode*);

  /// evaluate hits
  void evaluate_hits();

  /// evaluate clusters
  void evaluate_clusters();

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits(TrkrDefs::cluskey) const;

  /// true if histograms are initialized
  bool m_initialized = false;

  //! micromegas geometry
  PHG4CylinderGeomContainer* m_micromegas_geonode = nullptr;

  /// Acts tracking geometry for surface lookup
  ActsGeometry* m_tGeometry = nullptr;

  /// cluster map
  TrkrClusterContainer* m_cluster_map = nullptr;

  /// hitsets
  TrkrHitSetContainer* m_hitsets = nullptr;

  /// clusters to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

  /// hit to g4hit association
  TrkrHitTruthAssoc* m_hit_truth_map = nullptr;

  /// g4 hits
  PHG4HitContainer* m_g4hits_micromegas = nullptr;

  PHG4TruthInfoContainer* m_truthContainer;
  /// list of relevant layers
  /* it is filled at Init stage. It should not change for the full run */
  std::set<int> m_layers;
  ClusterErrorPara _ClusErrPara;
  int m_cluster_version = 4;
};

#endif
