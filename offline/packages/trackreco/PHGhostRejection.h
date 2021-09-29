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
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TpcSeedTrackMap;

class PHGhostRejection : public SubsysReco
{
 public:

  PHGhostRejection(const std::string &name = "PHGhostRejection");

  ~PHGhostRejection() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:

  int GetNodes(PHCompositeNode* topNode);
  void findGhostTracks();
  bool checkClusterSharing(SvtxTrack *tr1, SvtxTrack *tr2);

SvtxTrackMap *_track_map{nullptr};

  double _phi_cut = 0.01;
  double _eta_cut = 0.004;
  double _x_cut = 0.3;
  double _y_cut = 0.3;
  double _z_cut = 0.4;

};

#endif // PHGHOSTREJECTION_H
