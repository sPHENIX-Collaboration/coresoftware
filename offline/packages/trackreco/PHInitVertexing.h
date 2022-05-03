/*!
 *  \file		PHInitVertexing.h
 *  \brief		Base class for inital vertexing
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHINITVERTEXING_H
#define TRACKRECO_PHINITVERTEXING_H

// PHENIX includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <string>

// forward declarations
class PHCompositeNode;

class TrkrClusterContainer;
class SvtxVertexMap;

/// \class PHInitVertexing
///
/// \brief Base class for inital vertexing
///
class PHInitVertexing : public SubsysReco
{
 public:
  PHInitVertexing(const std::string &name = "PHInitVertexing");
  ~PHInitVertexing() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

 protected:

  /// setup interface for trackers, called in InitRun, setup things like pointers to nodes.
  /// overrided in derived classes
  virtual int Setup(PHCompositeNode *topNode);

  /// process event interface for trackers, called in process_event.
  /// implemented in derived classes
  virtual int Process(PHCompositeNode *topNode) = 0;


  TrkrClusterContainer *_cluster_map = nullptr;
  SvtxVertexMap *_vertex_map = nullptr;

 private:
  /// create new node output pointers
  int CreateNodes(PHCompositeNode *topNode);

  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);
};

#endif
