#ifndef TRACKRECO_PHTRUTHCLUSTERING_HP_H
#define TRACKRECO_PHTRUTHCLUSTERING_HP_H

#include <fun4all/SubsysReco.h>
#include <phool/PHObject.h>
#include <trackbase/TrkrDefs.h>

#include <map>
#include <set>
#include <string>
#include <vector>

class PHG4Hit;
class PHG4HitContainer;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;

class PHTruthClustering_hp : public SubsysReco
{
  public:

  /// constructor
  PHTruthClustering_hp( const std::string& = "PHTruthClustering_hp" );

  /// global initialization
  virtual int Init(PHCompositeNode*);

  /// run initialization
  virtual int InitRun(PHCompositeNode*);

  /// event processing
  virtual int process_event(PHCompositeNode*);

  /// end of processing
  virtual int End(PHCompositeNode*);

  private:

  /// load nodes
  int load_nodes( PHCompositeNode* );

  /// evaluate clusters
  void replace_clusters();

  // get geant hits associated to a cluster
  using G4HitSet = std::set<PHG4Hit*>;
  G4HitSet find_g4hits( TrkrDefs::cluskey ) const;

  // nodes
  TrkrClusterContainer* _cluster_map = nullptr;
  TrkrClusterHitAssoc* _cluster_hit_map = nullptr;
  TrkrHitTruthAssoc* _hit_truth_map = nullptr;

  PHG4HitContainer* _g4hits_tpc = nullptr;
  PHG4HitContainer* _g4hits_intt = nullptr;
  PHG4HitContainer* _g4hits_mvtx = nullptr;
  PHG4HitContainer* _g4hits_micromegas = nullptr;

};

#endif
