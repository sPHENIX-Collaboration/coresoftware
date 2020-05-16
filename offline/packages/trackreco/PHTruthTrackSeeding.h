/*!
 *  \file		PHTruthTrackSeeding.h
 *  \brief		Vertexing using truth info
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRUTHTRACKSEEDING_H
#define TRACKRECO_PHTRUTHTRACKSEEDING_H

#include "PHTrackSeeding.h"

#include <string>  // for string

// forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;

//class SvtxHitMap;
//class PHG4CellContainer;

/// \class PHTruthTrackSeeding
///
/// \brief Vertexing using truth info
///

class PHTruthTrackSeeding : public PHTrackSeeding
{
 public:
  PHTruthTrackSeeding(const std::string& name = "PHTruthTrackSeeding");
  virtual ~PHTruthTrackSeeding() {}

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
  int Setup(PHCompositeNode* topNode);

  int Process(PHCompositeNode* topNode);

  int End();

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode* topNode);

  PHG4TruthInfoContainer* _g4truth_container;

  PHG4HitContainer* phg4hits_tpc;
  PHG4HitContainer* phg4hits_intt;
  PHG4HitContainer* phg4hits_mvtx;

  TrkrHitTruthAssoc* hittruthassoc;
  TrkrClusterHitAssoc* clusterhitassoc;

  //SvtxHitMap* hitsmap;
  //PHG4CellContainer* cells_svtx;
  //PHG4CellContainer* cells_intt;
  //PHG4CellContainer* cells_maps;

  unsigned int _min_clusters_per_track;
  unsigned int _min_layer;
  unsigned int _max_layer;

  //! minimal truth momentum cut
  double _min_momentum;
};

#endif
