#ifndef HELIXHOUGH_VERTEXFINDER_H
#define HELIXHOUGH_VERTEXFINDER_H

#include <Eigen/Core>

// standard includes
#include <vector>

class SimpleTrack3D;

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
