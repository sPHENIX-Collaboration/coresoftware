#ifndef KALMAN_HELIXKALMANSTATE_H
#define KALMAN_HELIXKALMANSTATE_H

#include <Eigen/Core>

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
