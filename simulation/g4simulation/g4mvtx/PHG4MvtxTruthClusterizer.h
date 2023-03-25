#ifndef G4MVTX_TRUTHCLUSTERIZER__H
#define G4MVTX_TRUTHCLUSTERIZER__H

#include <g4tracking/TruthClusterizerBase.h>
  // This is the helper class for PHG4MvtxHitReco. It will cluster the PHG4Hits
  // which come from each individual embedded ("Truth") track. It uses pointers
  // to the modules (1) MvtxHitPruner, (2) PHG4MvtxDigitizer, and (3)
  // MvtxClusterizer, which also do the clustering for all the PHG4Hits
  // together (from all truth tracks together + all other tracks + noise).
  //
  // Because it makes calls to these other modules, t needs to be initialized
  // from the Fun4All driving macro *after* adding the above listed modules.

class PHG4MvtxDigitizer;
class MvtxHitPruner;
class MvtxClusterizer;
class PHCompositeNode;

class  PHG4MvtxTruthClusterizer : public TruthClusterizerBase {
  public:
  PHCompositeNode*   m_topNode     { nullptr };
  PHG4MvtxDigitizer* m_digitiser   ;
  MvtxHitPruner*     m_pruner      ;
  MvtxClusterizer*   m_clusterizer ;

  PHG4MvtxTruthClusterizer (
        PHG4MvtxDigitizer* _digitiser
      , MvtxHitPruner*     _pruner
      , MvtxClusterizer*   _clusterizer
      , int                _verbosity = 0 );

  void init_run(PHCompositeNode*& _topNode);
  int clusterize_hits(TrkrClusterContainer*);

  void check_g4hit(PHG4Hit*);
  void end_of_event();
  ~PHG4MvtxTruthClusterizer(){};

  /* fn_clusterer = mvtx_clusterize; */
  /* int check_g4hit (PHG4Hit* hit, int (*fn_clusterer)() ) { return _check_g4hit(hit, fn_clusterer=mvtx_clusterize); }; */
  /* void end_process_event () { mvtx_clusterize }; */

  /* void end_process_event() { return _end_process_event(mvtx_clusterize(m_topNode, m_digitiser, m_pruner, m_clusterizer, m_hits, m_verbosity, m_clusters_pass) ); }; */
};

#endif
