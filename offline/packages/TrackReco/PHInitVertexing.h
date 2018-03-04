/*!
 *  \file		PHInitVertexing.h
 *  \brief		Base class for inital vertexing
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __H_PHInitVertexing_H__
#define __H_PHInitVertexing_H__

// PHENIX includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <string>

// forward declarations
class PHCompositeNode;

class SvtxClusterMap;
class SvtxVertexMap;

/// \class PHInitVertexing
///
/// \brief Projects into calorimeters and fills track cal fields
///
class PHInitVertexing : public SubsysReco
{

 public:

  PHInitVertexing(const std::string &name = "PHInitVertexing");
  virtual ~PHInitVertexing() {}

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  virtual int Setup(PHCompositeNode *topNode);
  virtual int Process() = 0;

 protected:

  SvtxClusterMap *_cluster_map;
  SvtxVertexMap *_vertex_map;

 private:

	/// create new node output pointers
	int CreateNodes(PHCompositeNode *topNode);

	/// fetch node pointers
	int GetNodes(PHCompositeNode *topNode);
};

#endif // __H_PHInitVertexing_H__
