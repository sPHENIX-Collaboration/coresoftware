#ifndef __TRACK3D_H__
#define __TRACK3D_H__

#include <vector>
#include "Cluster3D.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

class Track3D {
 public:
  Track3D();
  Track3D(float _phi, float _d, float _kappa, float _dzdl, float _z0);
  virtual  ~Track3D() {}

  float fit_track(float scale = 1.);
  void set_vertex_id(unsigned int vtx_id){vertex_id = vtx_id;}
  void reset();

  std::vector<Cluster3D> hits;
  std::vector<unsigned int> cluster_ids;
  float phi, d, kappa, dzdl, z0;
  unsigned int index;
  unsigned int vertex_id;
};

#endif // __TRACK3D__
