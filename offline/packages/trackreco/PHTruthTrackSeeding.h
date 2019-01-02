/*!
 *  \file		PHTruthTrackSeeding.h
 *  \brief		Vertexing using truth info
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __H_PHTruthTrackSeeding_H__
#define __H_PHTruthTrackSeeding_H__

#include <set>
#include <vector>
#include "PHTrackSeeding.h"

// forward declarations
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class SvtxHitMap;
class PHG4CellContainer;

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

  const std::set<unsigned int>& get_seeding_layers() const
  {
    return _seeding_layers;
  }

  void set_seeding_layers(const unsigned int a[], const unsigned int n)
  {
    _seeding_layers.clear();
    for (unsigned int i = 0; i < n; ++i) _seeding_layers.insert(a[i]);
  }

 protected:
  int Setup(PHCompositeNode* topNode);

  int Process();

 private:
  /// create new node output pointers
  int CreateNodes(PHCompositeNode* topNode);

  /// fetch node pointers
  int GetNodes(PHCompositeNode* topNode);

  PHG4TruthInfoContainer* _g4truth_container;

  PHG4HitContainer* phg4hits_svtx;
  PHG4HitContainer* phg4hits_intt;
  PHG4HitContainer* phg4hits_maps;

  SvtxHitMap* hitsmap;

  PHG4CellContainer* cells_svtx;
  PHG4CellContainer* cells_intt;
  PHG4CellContainer* cells_maps;

  /// seeding layers
  std::set<unsigned int> _seeding_layers;

  unsigned int _min_clusters_per_track;
};

#endif  //__H_PHTruthTrackSeeding_H__
