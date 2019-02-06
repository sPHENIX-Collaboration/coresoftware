/*!
 *  \file		PHTruthVertexing.h
 *  \brief		Vertexing using truth info
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRUTHVERTEXING_H
#define TRACKRECO_PHTRUTHVERTEXING_H

#include "PHInitVertexing.h"

#include <vector>

// forward declarations
class PHG4TruthInfoContainer;

/// \class PHTruthVertexing
///
/// \brief Vertexing using truth info
///

class PHTruthVertexing : public PHInitVertexing
{
 public:
  PHTruthVertexing(const std::string &name = "PHTruthVertexing");
  virtual ~PHTruthVertexing() {}

  void set_vertex_error(const float &x_err, const float &y_err, const float &z_err)
  {
    _vertex_error.resize(3);
    _vertex_error[0] = x_err;
    _vertex_error[1] = y_err;
    _vertex_error[2] = z_err;
  }

  const std::vector<float> &get_vertex_error() const
  {
    return _vertex_error;
  }

 protected:
  int Setup(PHCompositeNode *topNode);

  int Process();

 private:
  /// create new node output pointers
  int CreateNodes(PHCompositeNode *topNode);

  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);

  PHG4TruthInfoContainer *_g4truth_container;

  /// manually assigned vertex error (standard dev), cm
  std::vector<float> _vertex_error;
};

#endif
