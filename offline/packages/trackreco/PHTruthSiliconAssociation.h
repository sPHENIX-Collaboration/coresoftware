// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHTRUTHSILICONASSOCIATION_H
#define PHTRUTHSILICONASSOCIATION_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <set>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class PHG4Particle;
class TpcSeedTrackMap;
class TrkrClusterCrossingAssoc;

class PHTruthSiliconAssociation : public SubsysReco
{
 public:

  PHTruthSiliconAssociation(const std::string &name = "PHTruthSiliconAssociation");

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;
  
 private:

  int GetNodes(PHCompositeNode* topNode);
  void copySiliconClustersToCorrectedMap( );

  std::vector<PHG4Particle*> getG4PrimaryParticle(SvtxTrack *track);
  std::set<TrkrDefs::cluskey> getSiliconClustersFromParticle(PHG4Particle* g4particle);
  std::vector<short int> getInttCrossings(SvtxTrack *si_track);

  PHG4TruthInfoContainer* _g4truth_container{nullptr};
  PHG4HitContainer *_g4hits_tpc{nullptr};
  PHG4HitContainer *_g4hits_mvtx{nullptr};
  PHG4HitContainer *_g4hits_intt{nullptr};
  
  TrkrClusterContainer *_cluster_map{nullptr};
  TrkrClusterContainer *_corrected_cluster_map{nullptr};
  TrkrClusterHitAssoc *_cluster_hit_map{nullptr};
  TrkrHitTruthAssoc *_hit_truth_map{nullptr};
  SvtxTrackMap *_track_map{nullptr};
  SvtxTrack *_tracklet{nullptr};
  SvtxVertexMap * _vertex_map{nullptr};
  TpcSeedTrackMap *_seed_track_map{nullptr};
 TrkrClusterCrossingAssoc *_cluster_crossing_map{nullptr};

 std::string _tpcseed_track_map_name = "TpcSeedTrackMap";

};

#endif // PHTRUTHSILICONASSOCIATION_H
