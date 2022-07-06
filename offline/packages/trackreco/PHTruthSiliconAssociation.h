// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHTRUTHSILICONASSOCIATION_H
#define PHTRUTHSILICONASSOCIATION_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <set>
#include <vector>
#include <memory>
#include <gsl/gsl_rng.h>

class PHCompositeNode;
class TrackSeedContainer;
class TrackSeed;
class SvtxVertexMap;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class PHG4Particle;
class TrkrClusterCrossingAssoc;
class ActsGeometry;

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
  //  void copySiliconClustersToCorrectedMap( );
  //void makeSvtxSeedMap();

  unsigned int buildTrackSeed(std::set<TrkrDefs::cluskey> clusters, PHG4Particle *g4particle, TrackSeedContainer* container);

  std::vector<PHG4Particle*> getG4PrimaryParticle(TrackSeed *track);
  std::set<TrkrDefs::cluskey> getSiliconClustersFromParticle(PHG4Particle* g4particle);
  std::set<short int> getInttCrossings(TrackSeed *si_track) const;

  PHG4TruthInfoContainer* _g4truth_container{nullptr};
  PHG4HitContainer *_g4hits_tpc{nullptr};
  PHG4HitContainer *_g4hits_mvtx{nullptr};
  PHG4HitContainer *_g4hits_intt{nullptr};
  
  TrkrClusterContainer *_cluster_map{nullptr};
  TrkrClusterContainer *_corrected_cluster_map{nullptr};
  TrkrClusterHitAssoc *_cluster_hit_map{nullptr};
  TrkrHitTruthAssoc *_hit_truth_map{nullptr};
  TrackSeedContainer *_tpc_track_map{nullptr};
  TrackSeedContainer *_silicon_track_map{nullptr};
  TrackSeedContainer *_svtx_seed_map{nullptr};
  TrackSeed *_tracklet{nullptr};
  SvtxVertexMap * _vertex_map{nullptr};
  TrkrClusterCrossingAssoc *_cluster_crossing_map{nullptr};
  ActsGeometry *_tgeometry{nullptr};

  //! rng de-allocator
  class Deleter
  {
    public:
    //! deletion operator
    void operator() (gsl_rng* rng) const { gsl_rng_free(rng); }
  };

 //! random generator that conform with sPHENIX standard
  /*! using a unique_ptr with custom Deleter ensures that the structure is properly freed when parent object is destroyed */
  std::unique_ptr<gsl_rng, Deleter> m_rng;

};

#endif // PHTRUTHSILICONASSOCIATION_H
