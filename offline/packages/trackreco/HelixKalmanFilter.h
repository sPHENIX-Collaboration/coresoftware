#ifndef G4HOUGH_HELIXKALMANFILTER_H
#define G4HOUGH_HELIXKALMANFILTER_H

#include <Eigen/Core>

#include <vector>

class SimpleHit3D;
class HelixTrackState;

class HelixKalmanFilter 
{
public:
  HelixKalmanFilter(std::vector<float>& detector_radii, 
		       std::vector<float>& detector_material, float B);
  virtual ~HelixKalmanFilter() {};

  void addHit(SimpleHit3D& hit, HelixTrackState& state);
    
protected:
  void calculateProjections(SimpleHit3D& hit, HelixTrackState& state, 
		Eigen::Matrix<float,2,5>& H, Eigen::Matrix<float,2,1>& ha);
  void calculateMeasurements(SimpleHit3D& hit, Eigen::Matrix<float,2,1>& m, 
		Eigen::Matrix<float,2,2>& G);
  void calculateMSCovariance(HelixTrackState& state, Eigen::Matrix<float,5,5>& Q);
  void subtractProjections(Eigen::Matrix<float,2,1>& m, Eigen::Matrix<float,2,1>& ha, 
		Eigen::Matrix<float,2,1>& diff);
  void updateIntersection(HelixTrackState& state, int layer);  

private:
  // Functions called by calculateMSCovariance
  bool calculateScatteringVariance(HelixTrackState& state, float& var);
  void calculate_dbdt(Eigen::Matrix<float,3,2>& dbdt_out);
  void calculate_dpdb(Eigen::Vector3f& p, Eigen::Matrix<float,3,3>& dpdb);
  void calculate_dApdp(HelixTrackState& state, Eigen::Matrix<float,3,3>& dApdp, 
		Eigen::Vector3f& p, float phi, float cosphi, float sinphi);
  void calculate_dAdAp(HelixTrackState& state, Eigen::Matrix<float,5,3>& dAdAp, 
		float& phi_p, float& cosphi_p, float& sinphi_p);
  // Functions called by calculateProjections
  void calculate_dxda(SimpleHit3D& hit, HelixTrackState& state, 
		Eigen::Matrix<float,3,5>& dxda, float& x, float& y, float& z);
  std::vector<float> det_radii;
  std::vector<float> det_scatter_variance;
  unsigned int nlayers;
  float signk_store;
  float Bfield, Bfield_inv;
    

};


#endif
