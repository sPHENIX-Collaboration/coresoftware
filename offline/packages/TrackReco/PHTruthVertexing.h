/*!
 *  \file		PHInitVertexing.h
 *  \brief		Base class for inital vertexing
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __H_PHTruthVertexing_H__
#define __H_PHTruthVertexing_H__

#include "PHInitVertexing.h"

#include <vector>

// forward declarations
class PHG4TruthInfoContainer;


/// \class PHInitVertexing
///
/// \brief Projects into calorimeters and fills track cal fields
///


class PHTruthVertexing : public PHInitVertexing {

public:

	PHTruthVertexing(const std::string &name = "PHTruthVertexing");
  virtual ~PHTruthVertexing() {}

	int Setup(PHCompositeNode *topNode);

	int Process();

	const std::vector<float>& get_vertex_error() const {
		return _vertex_error;
	}

	void set_vertex_error(const float & x_err, const float & y_err, const float & z_err) {
		_vertex_error.resize(3);
		_vertex_error[0] = x_err;
		_vertex_error[1] = y_err;
		_vertex_error[2] = z_err;
	}

private:
	/// create new node output pointers
	int CreateNodes(PHCompositeNode *topNode);

	/// fetch node pointers
	int GetNodes(PHCompositeNode *topNode);

	PHG4TruthInfoContainer* _g4truth_container;

	/// manually assigned vertex error (standard dev), cm
	std::vector<float> _vertex_error;

};

#endif //__H_PHTruthVertexing_H__
