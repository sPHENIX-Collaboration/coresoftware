/*!
 *  \file		  PH3DVertexing.h
 *  \brief		Base class for track seeding
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PH3DVERTEXING_H
#define TRACKRECO_PH3DVERTEXING_H

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

/// \class PH3DVertexing
///
/// \brief Base class for track seeding
///
class PH3DVertexing : public SubsysReco
{
 public:
  PH3DVertexing(const std::string &name = "PH3DVertexing");
  virtual ~PH3DVertexing() {}

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

#endif
