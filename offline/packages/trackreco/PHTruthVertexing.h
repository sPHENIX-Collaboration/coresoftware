/*!
 *  \file		PHTruthVertexing.h
 *  \brief		Vertexing using truth info
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef TRACKRECO_PHTRUTHVERTEXING_H
#define TRACKRECO_PHTRUTHVERTEXING_H


// rootcint barfs with this header so we need to hide it
#include <gsl/gsl_rng.h>

#include <string>             // for string
#include <vector>

#include "PHInitVertexing.h"

// forward declarations
class PHCompositeNode;
class PHG4TruthInfoContainer;

/// \class PHTruthVertexing
///
/// \brief Vertexing using truth info
///

class PHTruthVertexing : public PHInitVertexing
{
 public:
  PHTruthVertexing(const std::string &name = "PHTruthVertexing");
  ~PHTruthVertexing() override;

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
  void associate_tracks(bool associate_tracks)
  {
    _associate_tracks = associate_tracks;
  }
  void set_embed_only(bool embed_only)
  {
    _embed_only = embed_only;
  }
  void set_track_map_name(std::string& name)
  { _track_map_name = name; }
 protected:

  int Setup(PHCompositeNode *topNode) override;

  int Process(PHCompositeNode *topNode) override;

  int End(PHCompositeNode * /*topNode*/) override;

 private:
  /// fetch node pointers
  int GetNodes(PHCompositeNode *topNode);

  void assignTracksVertices(PHCompositeNode *topNode);

  PHG4TruthInfoContainer *_g4truth_container;

  /// manually assigned vertex error (standard dev), cm
  std::vector<float> _vertex_error;

  bool _embed_only;
  bool _associate_tracks = false;
  std::string _track_map_name = "SvtxTrackMap";
  gsl_rng *m_RandomGenerator;
};

#endif
