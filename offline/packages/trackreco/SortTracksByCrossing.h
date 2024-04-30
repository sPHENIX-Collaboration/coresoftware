// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 *  \file		  PHSimpleVertexFinder
 *  \brief		Class for determining primary vertices usin g fitted tracks
 *  \author	 Tony Frawley <afrawley@fsu.edu>
 */

#ifndef SORTTRACKSBYCROSSING_H
#define SORTTRACKSBYCROSSING_H

#include <trackbase/TrackVertexCrossingAssoc_v1.h>

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <map>
#include <set>

#include <Eigen/Dense>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class TrkrCluster;
class TrackVertexCrossingAssoc;

class SortTracksByCrossing : public SubsysReco
{
 public:

  SortTracksByCrossing(const std::string &name = "SortTracksByCrossing");

  ~SortTracksByCrossing() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:

  int GetNodes(PHCompositeNode* topNode);
  int CreateNodes(PHCompositeNode* topNode);
  
  SvtxTrackMap *_track_map{nullptr};
  
  TrackVertexCrossingAssoc* _track_vertex_crossing_map;

  
};

#endif // SORTTRACKSBYCROSSING_H
