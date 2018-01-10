#include "HelixTrackState.h"

using namespace Eigen;


HelixTrackState::HelixTrackState() : position(0), phi(0.), d(0.), kappa(0.), z0(0.), dzdl(0.), nu(0.), x_int(0.), y_int(0.), z_int(0.), chi2(0.), C(Matrix<float,5,5>::Identity(5,5))
{
  
}


HelixTrackState::~HelixTrackState()
{
  
}
