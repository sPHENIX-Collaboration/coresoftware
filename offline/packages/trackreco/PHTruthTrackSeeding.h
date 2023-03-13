/*!
 *  \file		PHTruthTrackSeeding.h
 *  \brief		Vertexing using truth info
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRUTHTRACKSEEDING_H
#define TRACKRECO_PHTRUTHTRACKSEEDING_H

#include "PHTrackSeeding.h"
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <string>  // for string
#include <vector>
#include <memory>
#include <gsl/gsl_rng.h>

// forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class TrkrHitTruthAssoc;
class TrkrClusterContainer;
class TrkrClusterCrossingAssoc;
class SvtxClusterEval;
class TrackSeed;
class TrackSeedContainer;
class PHG4Particle;

/// \class PHTruthTrackSeeding
///
/// \brief Vertexing using truth info
///

class PHTruthTrackSeeding : public PHTrackSeeding
{
 public:
  PHTruthTrackSeeding(const std::string& name = "PHTruthTrackSeeding");

  unsigned int get_min_clusters_per_track() const
  {
    return _min_clusters_per_track;
  }

  void set_min_clusters_per_track(unsigned int minClustersPerTrack)
  {
    _min_clusters_per_track = minClustersPerTrack;
  }

  void set_min_layer(unsigned int minLayer)
  {
    _min_layer = minLayer;
  }

  void set_max_layer(unsigned int maxLayer)
  {
    _max_layer = maxLayer;
  }

  //! minimal truth momentum cut
  double get_min_momentum() const
  {
    return _min_momentum;
  }

  //! minimal truth momentum cut
  void set_min_momentum(double m)
  {
    _min_momentum = m;
  }

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process(PHCompositeNode* topNode) override;

  int End() override;

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode* topNode);
  int CreateNodes(PHCompositeNode* topNode);

  void buildTrackSeed(std::vector<TrkrDefs::cluskey> clusters, 
		      PHG4Particle *g4particle, TrackSeedContainer* container);
  PHG4TruthInfoContainer* m_g4truth_container = nullptr;

  /// get crossing id from intt clusters associated to track
  /* this is a copy of the code in PHTruthSiliconAssociation */
  std::set<short int> getInttCrossings(TrackSeed*) const;

  TrkrClusterContainer *m_clusterMap = nullptr;
  TrkrClusterCrossingAssoc *m_cluster_crossing_map = nullptr;
  PHG4HitContainer* phg4hits_tpc = nullptr;
  PHG4HitContainer* phg4hits_intt = nullptr;
  PHG4HitContainer* phg4hits_mvtx = nullptr;
  PHG4HitContainer* phg4hits_micromegas = nullptr;

  TrkrHitTruthAssoc* hittruthassoc = nullptr;
  SvtxClusterEval* _clustereval;

  unsigned int _min_clusters_per_track = 3;
  unsigned int _min_layer = 0;
  unsigned int _max_layer = 60;

  //! minimal truth momentum cut (GeV)
  double _min_momentum = 50e-3;

  TrackSeedContainer *_track_map_silicon = nullptr;
  TrackSeedContainer *_track_map_combined = nullptr;

  ActsGeometry *tgeometry = nullptr;

  bool _circle_fit_seed = false;

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

#endif
