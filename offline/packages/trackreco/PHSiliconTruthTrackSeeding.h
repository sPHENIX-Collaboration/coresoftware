/*!
 *  \file		PHSiliconTruthTrackSeeding.h
 *  \brief		Vertexing using truth info
 *  \author		Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_PHSILICONTRUTHTRACKSEEDING_H
#define TRACKRECO_PHSILICONTRUTHTRACKSEEDING_H

#include "PHTrackSeeding.h"

#include <string>  

// forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4HitContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;


/// \class PHSiliconTruthTrackSeeding
///
/// \brief Vertexing using truth info
///

class PHSiliconTruthTrackSeeding : public PHTrackSeeding
{
 public:
  PHSiliconTruthTrackSeeding(const std::string& name = "PHSiliconTruthTrackSeeding");

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

  PHG4TruthInfoContainer* _g4truth_container = nullptr;

  PHG4HitContainer* phg4hits_tpc = nullptr;
  PHG4HitContainer* phg4hits_intt = nullptr;
  PHG4HitContainer* phg4hits_mvtx = nullptr;
  PHG4HitContainer* phg4hits_micromegas = nullptr;

  TrkrHitTruthAssoc* hittruthassoc = nullptr;
  TrkrClusterHitAssoc* clusterhitassoc = nullptr;

  unsigned int _min_clusters_per_track = 2;
  unsigned int _min_layer = 0;
  unsigned int _max_layer = 6;

  //! minimal truth momentum cut (GeV)
  double _min_momentum = 50e-3;
};

#endif
