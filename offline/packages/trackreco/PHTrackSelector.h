// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  PHTrackCleaner
 *  \brief		Class for deciding which track based on a given TPC seed is the best one
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef PHTRACKSELECTOR_H
#define PHTRACKSELECTOR_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <map>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrkrClusterIterationMapv1;

class PHTrackSelector : public SubsysReco
{
 public:

  PHTrackSelector(const std::string &name = "PHTrackSelector");

  ~PHTrackSelector() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void SetMinTPCClusters(int nhits) {min_tpc_clusters = nhits;}
  void SetMinINTTHits(int nhits) {min_intt_hits = nhits;}
  void SetMinMVTXHits(int nhits) {min_mvtx_hits = nhits;}
  void SetChi2NDFHits(float max) {max_chi2_ndf = max;}
  void SetIteration(int iter) {_n_iter = iter;}

  void SetTrackMapName(const std::string &map_name) { _track_map_name = map_name; }

 private:

  int GetNodes(PHCompositeNode* topNode);

  std::string _track_map_name;
  SvtxTrackMap *_track_map{nullptr};
  SvtxTrack *_track{nullptr};
  TrkrClusterIterationMapv1* _iteration_map = nullptr;

  int _n_iter = 1;
  unsigned int min_tpc_clusters = 35;
  unsigned int min_mvtx_hits = 2;
  unsigned int min_intt_hits = 1;
  float max_chi2_ndf = 30;

};

#endif // PHTRACKSELECTOR_H
