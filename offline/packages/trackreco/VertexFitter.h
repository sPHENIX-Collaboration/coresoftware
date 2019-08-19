#ifndef G4HOUGH_VERTEXFITTER_H
#define G4HOUGH_VERTEXFITTER_H

#include <Eigen/Core>

// standard includes
#include <vector>

class SimpleTrack3D;

/// \class VertexFitter
///
/// \brief A class to determine the vertex
///
/// This class incorporates Newton's method for gradient descent
/// against an expo-dca^2 function.
///
class VertexFitter {
 public:
  VertexFitter();
  virtual ~VertexFitter() {}

  bool findVertex(std::vector<SimpleTrack3D>& tracks,
                  std::vector<Eigen::Matrix<float, 5, 5> >& covariances,
                  std::vector<float>& vertex, float sigma, bool fix_xy = false);

  bool findVertex(std::vector<SimpleTrack3D>& tracks,
                  std::vector<float>& vertex, float sigma, bool fix_xy = false);

 protected:
};

#endif  // __VERTEXFITTER_H__
