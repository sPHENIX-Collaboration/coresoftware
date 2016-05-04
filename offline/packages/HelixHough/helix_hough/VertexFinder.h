#ifndef __VERTEXFINDER_H__
#define __VERTEXFINDER_H__

// Helix Hough includes
#include <SimpleTrack3D.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

// standard includes
#include <vector>

/// \class VertexFinder
///
/// \brief A class to determine the vertex
///
/// This class incorporates Newton's method for gradient descent
/// against an expo-dca^2 function.
///
class VertexFinder {
 public:
  VertexFinder();
  virtual ~VertexFinder() {}

  bool findVertex(std::vector<SimpleTrack3D>& tracks,
                  std::vector<Eigen::Matrix<float, 5, 5> >& covariances,
                  std::vector<float>& vertex, float sigma, bool fix_xy = false);

  bool findVertex(std::vector<SimpleTrack3D>& tracks,
                  std::vector<float>& vertex, float sigma, bool fix_xy = false);

 protected:
};

#endif  // __VERTEXFINDER_H__
