// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  PHGhostRejection
 *  \brief		Class for deciding which track based on a given TPC seed is the best one
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef PHGHOSTREJECTION_H
#define PHGHOSTREJECTION_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>

#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class TrackSeedContainer;
class TrkrCluster;
class TrackSeed;
class TrkrClusterContainer;

class PHGhostRejection
{
 public:
  PHGhostRejection() {}
  PHGhostRejection(unsigned int verbosity) : m_verbosity { verbosity } {};

  ~PHGhostRejection();

  void rejectGhostTracks(std::vector<float> &trackChi2, TrackSeedContainer* tracks, std::map<TrkrDefs::cluskey, Acts::Vector3>& positions);
  void verbosity(int verb) { m_verbosity = verb; }
  /* void trackSeedContainer(TrackSeedContainer *seeds) { m_trackMap = seeds; } */
  /* void positionMap(std::map<TrkrDefs::cluskey, Acts::Vector3> &map) { m_positions = map; } */

 private:
  bool checkClusterSharing(TrackSeed *tr1, TrackSeed *tr2);

  double _phi_cut = 0.01;
  double _eta_cut = 0.004;
  double _x_cut = 0.3;
  double _y_cut = 0.3;
  double _z_cut = 0.4;

  // cuts on minimally interesting tracks -- 
  // here to remove noise
  double _min_pt_cut = 0.2;
  bool   _must_span_sectors = true;
  size_t _min_clusters = 6;

  unsigned int m_verbosity = 0;

  /* TrackSeedContainer *m_trackMap = nullptr; */

  /* std::map<TrkrDefs::cluskey, Acts::Vector3> m_positions; */
};

#endif  // PHGHOSTREJECTION_H
