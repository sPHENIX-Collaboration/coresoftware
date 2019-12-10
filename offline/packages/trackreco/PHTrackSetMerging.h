/*!
 *  \file		PHTrackSetMerging.h
 *  \brief		Base class for track container merging
 *  \author		Christof Roland <cer@mit.edu>
 */

#ifndef TRACKRECO_PHTRACKSETMERGING_H
#define TRACKRECO_PHTRACKSETMERGING_H

// PHENIX includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <string>

// forward declarations
class PHCompositeNode;

//class SvtxClusterMap;
class TrkrClusterContainer;
class SvtxVertexMap;
class SvtxTrackMap;
class AssocInfoContainer;

/// \class PHTrackSetMerging
///
/// \brief Base class for track seeding
///
class PHTrackSetMerging : public SubsysReco
{
 public:
  PHTrackSetMerging(const std::string &name = "PHTrackSetMerging");
  virtual ~PHTrackSetMerging() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_track_map_name_in1(const std::string &map_name) { _track_map_name_in1 = map_name; }
  void set_track_map_name_in2(const std::string &map_name) { _track_map_name_in2 = map_name; }
  void set_track_map_name_out(const std::string &map_name) { _track_map_name_out = map_name; }

 protected:
  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// overrided in derived classes
  virtual int Setup(PHCompositeNode *topNode);

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  virtual int Process() = 0;

  /// Called in SubsysReco::End
  virtual int End() = 0;

  //SvtxClusterMap *_cluster_map;
  TrkrClusterContainer *_cluster_map;
  SvtxVertexMap *_vertex_map;
  SvtxTrackMap *_track_map_in1;
  SvtxTrackMap *_track_map_in2;
  SvtxTrackMap *_track_map_out;
  AssocInfoContainer *_assoc_container;

  std::string _track_map_name_in1;
  std::string _track_map_name_in2;
  std::string _track_map_name_out;

 private:
  /// create new node output pointers
  int CreateNodes(PHCompositeNode *topNode);

  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);
};

#endif
