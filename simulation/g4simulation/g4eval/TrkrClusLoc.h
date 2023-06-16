#ifndef TRKRCLUSLOC__H
#define TRKRCLUSLOC__H
// A few simple structures that are used in the evaluation and matching of clusters
// Plus a few functions for working with clusters

#include <trackbase/TrkrDefs.h>
#include <Eigen/Core>

struct TrkrClusLoc {
  int layer{0};
  Eigen::Vector3d gloc{};
  float phi{0};
  float phisize{0};
  float z{0};
  float zsize{0};
  TrkrDefs::cluskey ckey{};
};

#endif
