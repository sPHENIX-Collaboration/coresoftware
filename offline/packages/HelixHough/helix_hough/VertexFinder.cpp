
#include "VertexFinder.h"

// Helix Hough includes
#include <VertexFitFunc.h>

// FitNewton includes
#include <NewtonMinimizerGradHessian.h>

// Eigen includes
#include <Eigen/Core>

class SimpleTrack3D;

using namespace std;

using namespace FitNewton;
using namespace Eigen;

VertexFinder::VertexFinder()
{

}

/// The method will find the vertex given an initial guess and a list of
/// track candidates.
///
/// \param[in] tracks List of track candidates
/// \param[in,out] vertex The vertex position in x,y,z
///
/// \return true if successful
///
bool VertexFinder::findVertex(vector<SimpleTrack3D>& tracks, vector<float>& vertex, float sigma, bool fix_xy)
{
  vector<Matrix<float,5,5> > covariances;
  return findVertex(tracks, covariances, vertex, sigma, fix_xy);
}



bool VertexFinder::findVertex(vector<SimpleTrack3D>& tracks, vector<Matrix<float,5,5> >& covariances, vector<float>& vertex, float sigma, bool fix_xy)
{
  VertexFitFunc _vertexfit;
  FitNewton::NewtonMinimizerGradHessian _minimizer;

  // setup function to minimize
  // expo-dca2 => ~dca^2 w/ extreme outlier de-weighting
  _vertexfit.setTracks(&tracks);
  _vertexfit.setCovariances(&covariances);
    
  // setup the minimizer
  // fast Newton gradient minimization
  _minimizer.setFunction(&_vertexfit);
    
  // set an initial guess and uncertainties
  VectorXd start_point = VectorXd::Zero(vertex.size()); // input initial guess
  for(unsigned int i=0;i<vertex.size();++i)
  {
    start_point(i) = vertex[i];
  }
  if(fix_xy)
  {
    _minimizer.fixParameter(0);
    _minimizer.fixParameter(1);
  }
  _vertexfit.setFixedPar(0, sigma);                     // index = 0 used for x,y,z uncertainties

  // output storage
  VectorXd min_point = VectorXd::Zero(vertex.size());   // output vector for minimize method below

  // minimize
  _minimizer.minimize(start_point, min_point, 1.0e-12, 48, 1.0e-18);

  // store output vertex spatial point
  for(unsigned i=0; i<vertex.size(); i++)
    {
      vertex[i] = min_point(i);
    }

  return true;
}
