#ifndef __SIMPLETRACK3D__
#define __SIMPLETRACK3D__

#include <vector>
#include "SimpleHit3D.h"

class SimpleTrack3D {
 public:
  SimpleTrack3D();
  ~SimpleTrack3D() {}

  float fit_track(float scale = 1.);
  void set_vertex_id(unsigned int vtx_id){vertex_id = vtx_id;}
  void reset();

  std::vector<SimpleHit3D> hits;
  std::vector<unsigned int> cluster_ids;
  float phi, d, kappa, dzdl, z0;
  unsigned int index;
  unsigned int vertex_id;
};

#endif // __SIMPLETRACK3D__
