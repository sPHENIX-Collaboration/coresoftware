#ifndef KALMAN_HELIXKALMAN_H
#define KALMAN_HELIXKALMAN_H

#include <Eigen/Core>

class HelixKalmanState;
class SimpleHit3D;


class HelixKalman
{
  public:
    HelixKalman(float B);
    virtual ~HelixKalman();
    
    void addHit(SimpleHit3D& hit, HelixKalmanState& state);
    
  protected:
    virtual void calculateProjections(SimpleHit3D& hit, HelixKalmanState& state, Eigen::Matrix<float,2,5>& H, Eigen::Matrix<float,2,1>& ha) = 0;
    virtual void calculateMeasurements(SimpleHit3D& hit, Eigen::Matrix<float,2,1>& m, Eigen::Matrix<float,2,2>& G) = 0;
    virtual bool calculateScatteringVariance(HelixKalmanState& state, float& var) = 0;
    virtual void updateIntersection(HelixKalmanState& state, int layer) = 0;
    virtual void subtractProjections(Eigen::Matrix<float,2,1>& m, Eigen::Matrix<float,2,1>& ha, Eigen::Matrix<float,2,1>& diff);
    
    void calculateMSCovariance(HelixKalmanState& state, Eigen::Matrix<float,5,5>& Q);
    void calculate_dbdt(Eigen::Matrix<float,3,2>& dbdt_out);
    void calculate_dpdb(Eigen::Vector3f& p, Eigen::Matrix<float,3,3>& dpdb);
    void calculate_dApdp(HelixKalmanState& state, Eigen::Matrix<float,3,3>& dApdp, Eigen::Vector3f& p, float phi, float cosphi, float sinphi);
    void calculate_dAdAp(HelixKalmanState& state, Eigen::Matrix<float,5,3>& dAdAp, float& phi_p, float& cosphi_p, float& sinphi_p);
    
    float Bfield, Bfield_inv;
};


#endif
