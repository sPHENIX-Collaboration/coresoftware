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
class TpcSeedTrackMap;
class TrackSeed;

class PHGhostRejection 
{
 public:

  PHGhostRejection(unsigned int verbosity );

  ~PHGhostRejection();

  void rejectGhostTracks(TrackSeedContainer *trackMap, std::vector<float>& trackChi2);
  void setVerbosity(int verb) { m_verbosity = verb; }

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

};

#endif // PHGHOSTREJECTION_H
