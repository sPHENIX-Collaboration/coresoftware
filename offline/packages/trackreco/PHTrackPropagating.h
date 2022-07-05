/*!
 *  \file		  PHTrackPropagating.h
 *  \brief		Base class for track seeding
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRACKPROPAGATING_H
#define TRACKRECO_PHTRACKPROPAGATING_H

// PHENIX includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <string>

// forward declarations
class PHCompositeNode;

class TrkrClusterContainer;
class SvtxVertexMap;
class SvtxTrackMap;

/// \class PHTrackPropagating
///
/// \brief Base class for track seeding
///
class PHTrackPropagating : public SubsysReco
{
 public:
  PHTrackPropagating(const std::string &name = "PHTrackPropagating");
  ~PHTrackPropagating() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void SetUseTruthClusters(bool setit){_use_truth_clusters = setit;}

 protected:
  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// overrided in derived classes
  virtual int Setup(PHCompositeNode *topNode);

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  virtual int Process() = 0;

  ///
  virtual int End() = 0;


  //SvtxClusterMap *_cluster_map;
  TrkrClusterContainer *_cluster_map = nullptr;
  SvtxVertexMap *_vertex_map = nullptr;
  SvtxTrackMap *_track_map = nullptr;

  std::string _track_map_name = "SvtxTrackMap";

  bool _use_truth_clusters = false;

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);
};

#endif
