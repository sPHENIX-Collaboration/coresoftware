/*!
 *  \file		  PHTrackFitting.h
 *  \brief		Base class for track seeding
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRACKFITTING_H
#define TRACKRECO_PHTRACKFITTING_H

// PHENIX includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <string>

// forward declarations
class PHCompositeNode;

//class SvtxClusterMap;
class TrkrClusterContainer;
class TrkrHitSetContainer;
class SvtxVertexMap;
class SvtxTrackMap;

/// \class PHTrackFitting
///
/// \brief Base class for track seeding
///
class PHTrackFitting : public SubsysReco
{
 public:
  PHTrackFitting(const std::string &name = "PHTrackFitting");
  ~PHTrackFitting() override {}

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  //virtual const std::set<unsigned int> &get_seeding_layers() const = 0;

  //virtual void set_seeding_layers(const unsigned int a[], const unsigned int n) = 0;

  //void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }

 protected:
  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// overrided in derived classes
  virtual int Setup(PHCompositeNode *topNode);

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  virtual int Process() = 0;

  //SvtxClusterMap *_cluster_map;
  TrkrClusterContainer *_cluster_map;
  TrkrHitSetContainer  *_hitsets = nullptr;
  SvtxVertexMap *_vertex_map;
  SvtxTrackMap *_track_map;

  std::string _track_map_name;

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);
};

#endif
