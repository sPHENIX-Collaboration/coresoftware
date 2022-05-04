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

#include <string>
#include <vector>
#include <map>

class PHCompositeNode;
class TrackSeedContainer;
class TrkrCluster;
class TrackSeed;
class ActsSurfaceMaps;
class ActsTrackingGeometry;
class TrkrClusterContainer;

class PHGhostRejection 
{
 public:
  PHGhostRejection() {}
  PHGhostRejection(unsigned int verbosity);

  ~PHGhostRejection();

  void rejectGhostTracks(std::vector<float>& trackChi2);
  void verbosity(int verb) { m_verbosity = verb; }
  void trackSeedContainer(TrackSeedContainer *seeds){m_trackMap = seeds;}
  void clusterContainer(TrkrClusterContainer *clus) {m_clusters = clus;}
  void surfMaps(ActsSurfaceMaps *maps) {m_surfmaps = maps;}
  void geometry(ActsTrackingGeometry *geom) {m_tGeometry = geom;}

 private:

  bool checkClusterSharing(TrackSeed *tr1, unsigned int trid1,
			   TrackSeed *tr2, unsigned int trid2);

  double _phi_cut = 0.01;
  double _eta_cut = 0.004;
  double _x_cut = 0.3;
  double _y_cut = 0.3;
  double _z_cut = 0.4;
  int _n_iteration = 0;
  unsigned int m_verbosity = 0;
  
  TrackSeedContainer *m_trackMap = nullptr;
  TrkrClusterContainer *m_clusters = nullptr;
  ActsSurfaceMaps *m_surfmaps = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;


};

#endif // PHGHOSTREJECTION_H
