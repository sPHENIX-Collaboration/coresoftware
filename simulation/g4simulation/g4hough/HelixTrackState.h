#ifndef G4HOUGH_HELIXTRACKSTATE_H
#define G4HOUGH_HELIXTRACKSTATE_H

#include <Eigen/Core>

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
