#ifndef __HELIXKALMANSTATE__
#define __HELIXKALMANSTATE__

#include <Eigen/LU>
#include <Eigen/Core>
#include <Eigen/Dense>


class HelixKalmanState
{
  public:
    HelixKalmanState();
    ~HelixKalmanState();
    
    unsigned int position;
    float phi, d, kappa, z0, dzdl;
    // nu^2 = kappa
    float nu;
    float x_int, y_int, z_int;
    float chi2;
    Eigen::Matrix<float,5,5> C;
};

#endif
