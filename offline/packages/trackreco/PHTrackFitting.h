/*!
 *  \file		  PHTrackFitting.h
 *  \brief		Base class for track seeding
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __H_PHTrackFitting_H__
#define __H_PHTrackFitting_H__

// PHENIX includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <set>
#include <string>

// forward declarations
class PHCompositeNode;

class SvtxClusterMap;
class SvtxVertexMap;
class SvtxTrackMap;
class AssocInfoContainer;

/// \class PHTrackFitting
///
/// \brief Base class for track seeding
///
class PHTrackFitting : public SubsysReco
{
 public:
  PHTrackFitting(const std::string &name = "PHTrackFitting");
  virtual ~PHTrackFitting() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  virtual const std::set<unsigned int> &get_seeding_layers() const = 0;

  virtual void set_seeding_layers(const unsigned int a[], const unsigned int n) = 0;

 protected:
  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// overrided in derived classes
  virtual int Setup(PHCompositeNode *topNode);

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  virtual int Process() = 0;

  SvtxClusterMap *_cluster_map;
  SvtxVertexMap *_vertex_map;
  SvtxTrackMap *_track_map;
  AssocInfoContainer *_assoc_container;

 private:
  /// create new node output pointers
  int CreateNodes(PHCompositeNode *topNode);

  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);
};

#endif  // __H_PHTrackFitting_H__
