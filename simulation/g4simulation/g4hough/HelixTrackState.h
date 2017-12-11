#ifndef __HELIXTRACKSTATE_H__
#define __HELIXTRACKSTATE_H__

#include <Eigen/LU>
#include <Eigen/Core>
#include <Eigen/Dense>


class HelixTrackState
{
  public:
    HelixTrackState();
    ~HelixTrackState();
    
    unsigned int position;
    float phi, d, kappa, z0, dzdl;
    // nu^2 = kappa
    float nu;
    float x_int, y_int, z_int;
    float chi2;
    Eigen::Matrix<float,5,5> C;
};

#endif
