// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  PHTrackCleaner
 *  \brief		Class for deciding which track based on a given TPC seed is the best one
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef PHTRACKCLEANER_H
#define PHTRACKCLEANER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <map>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrackSeedContainer;

class PHTrackCleaner : public SubsysReco
{
 public:

  PHTrackCleaner(const std::string &name = "PHTrackCleaner");

  ~PHTrackCleaner() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 void set_pp_mode(const bool flag){_pp_mode = flag ;}

 private:

  int GetNodes(PHCompositeNode* topNode);
  void findGhostTracks();

SvtxTrackMap *_track_map{nullptr};
SvtxTrack *_track{nullptr};
 TrackSeedContainer *_tpc_seed_map{nullptr};

 double min_ndf = 25;
 bool _pp_mode = false;

};

#endif // PHTRACKCLEANER_H
