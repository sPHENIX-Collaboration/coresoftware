/*!
 *  \file		  PHTrackSeeding.h
 *  \brief		Base class for track seeding
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRACKSEEDING_H
#define TRACKRECO_PHTRACKSEEDING_H

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

/// \class PHTrackSeeding
///
/// \brief Base class for track seeding
///
class PHTrackSeeding : public SubsysReco
{
 public:
  PHTrackSeeding(const std::string &name = "PHTrackSeeding");
  virtual ~PHTrackSeeding() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }

  //virtual const std::set<unsigned int>& get_seeding_layers() const = 0;

  //virtual void set_seeding_layers(const unsigned int a[], const unsigned int n) = 0;

 protected:
  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// overrided in derived classes
  virtual int Setup(PHCompositeNode *topNode);

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  virtual int Process(PHCompositeNode *topNode) = 0;

  /// Called in SubsysReco::End
  virtual int End() = 0;

  //SvtxClusterMap *_cluster_map;
  TrkrClusterContainer *_cluster_map;
  SvtxVertexMap *_vertex_map;
  SvtxTrackMap *_track_map;
  AssocInfoContainer *_assoc_container;

  std::string _track_map_name;

 private:
  /// create new node output pointers
  int CreateNodes(PHCompositeNode *topNode);

  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);
};

#endif
