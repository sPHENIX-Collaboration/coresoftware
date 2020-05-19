#include "ActsCovarianceRotater.h"


Acts::BoundSymMatrix ActsCovarianceRotater::rotateSvtxTrackCovToActs(
	             const SvtxTrack *track)
{
  Acts::BoundSymMatrix matrix = Acts::BoundSymMatrix::Zero();

  const double x = track->get_x();
  const double y = track->get_y();

  const double px = track->get_px();
  const double py = track->get_py();
  const double pz = track->get_pz();
  const double p = sqrt(px * px + py * py + pz * pz);

  const double phiPos = atan2(y,x);
  const int charge = track->get_charge();

  Acts::BoundSymMatrix svtxTrackCov = Acts::BoundSymMatrix::Zero();
  
  for (int i = 0; i < 6; ++i)
    {
      for (int j = 0; j < 6; ++j)
	{
	  svtxTrackCov(i,j) = track->get_error(i,j);
	}
    }

  printMatrix("Initial covariance :", svtxTrackCov);

  /// Construct the jacobian transformation matrix. Can be determined
  /// by writing Acts quantities in terms of (x,y,z,px,py,pz) and taking
  /// derivatives
  /// Rotation matrix does not need units because everything is defined
  /// in a unitless way, intentionally
  Acts::BoundSymMatrix rotation = Acts::BoundSymMatrix::Zero();

  /// Make a unit p vector for the rotation
  const double uPx = px / p;
  const double uPy = py / p;
  const double uPz = pz / p;
  const double uP = sqrt(uPx * uPx + uPy * uPy + uPz * uPz);

  /// d0 = sqrt(x*x + y*y) when reference point is at (0,0,0)
  rotation(0,0) = cos(phiPos);
  rotation(0,1) = sin(phiPos);

  /// z to z0 rotation is trivial
  rotation(1,2) = 1;
  
  /// Rotating px,py,pz to phi, theta, p is just a rotation from cartesian to spherical
  /// Rotate to phi
  rotation(2,3) = -1 * uPy / (uPx * uPx + uPy * uPy);
  rotation(2,4) = uPx / (uPx * uPx + uPy * uPy);

  /// Rotate to theta. Leave uP in for clarity even though it is trivially 1
  rotation(3,3) = uPx * uPz / (uP* uP * sqrt(uPx * uPx + uPy * uPy));
  rotation(3,4) = uPy * uPz / (uP * uP * sqrt(uPx * uPx + uPy * uPy));
  rotation(3,5) = -1 * sqrt(uPx * uPx + uPy * uPy) / (uP * uP);

  /// Rotate to p. Since p is a unit vector it is trivially 1 but leave it
  /// in for clarity
  rotation(4,3) = uPx / uP;
  rotation(4,4) = uPy / uP;
  rotation(4,5) = uPz / uP;

  /// time component is 0, so we have no entries for rotation(5,j)

  printMatrix("Rotation matrix is : ", rotation);

  matrix = rotation * svtxTrackCov * rotation.transpose();

  printMatrix("Rotated matrix is : ", matrix);

  /// Now rotate to q/p, which is mostly trivial
  Acts::BoundSymMatrix qprotation = Acts::BoundSymMatrix::Zero();
  qprotation(0,0) = 1;
  qprotation(1,1) = 1;
  qprotation(2,2) = 1;
  qprotation(3,3) = 1;
  /// var(q/p) = (d(1/p)/dp)^2 * var(p) = (-1/p^2)^2 * var(p), from acts devel
  qprotation(4,4) = charge * charge / (p * p * p * p);
  qprotation(5,5) = 1;
  
  matrix = qprotation * matrix * rotation.transpose();

  /// Need to make sure units are what acts expects. Acts defaults to mm and GeV
  /// so covariance entries with position component must be multiplied by 10
  for(int i = 0; i < 6; ++i)
    {
      for(int j = 0; j < 6; ++j)
	{
	  if(i < 2 && j < 2)
	    matrix(i,j) *= Acts::UnitConstants::cm2;
	  else if((i < 2 && j < 5) || (i < 5 && j < 2))
	    matrix(i,j) *= Acts::UnitConstants::cm;	  
	}
    }

  return matrix;
}


Acts::BoundSymMatrix ActsCovarianceRotater::rotateActsCovToSvtxTrack(
                     const Acts::KalmanFitterResult<SourceLink>& fitOutput)
{

  Acts::BoundSymMatrix matrix = Acts::BoundSymMatrix::Zero();
  
  const auto& params = fitOutput.fittedParameters.value();
  auto covarianceMatrix = *params.covariance();
  
  printMatrix("Initial Acts covariance: ", covarianceMatrix);

  const double px = params.momentum()(0);
  const double py = params.momentum()(1);
  const double pz = params.momentum()(2);
  const double p = sqrt(px * px + py * py + pz * pz);
  
  const double x = params.position()(0);
  const double y = params.position()(1);

  const int charge = params.charge();
  const double phiPos = atan2(x, y);

  /// Return to sPHENIX units of cm
  for(int i = 0; i < 6; ++i)
    {
      for(int j = 0; j < 6; ++j)
	{
	  if(i < 2 && j < 2)
	    matrix(i,j) /= Acts::UnitConstants::cm2;
	  else if((i < 2 && j < 5) || (i < 5 && j < 2))
	    matrix(i,j) /= Acts::UnitConstants::cm;	  
	}
    }

  /// First rotate the covariance matrix from (d0,z0,phi,theta,q/p,t) 
  /// to (d0,z0,phi,theta,p,t) since it is easier to work with
  Acts::BoundSymMatrix qprotation = Acts::BoundSymMatrix::Zero();
  qprotation(0,0) = 1;
  qprotation(1,1) = 1;
  qprotation(2,2) = 1;
  qprotation(3,3) = 1;
  /// var(q/p) = (d(1/p)/dp)^2 * var(p) = (-1/p^2)^2 * var(p)
  qprotation(4,4) = (p * p * p * p) / (charge * charge);
  qprotation(5,5) = 1;

  covarianceMatrix = qprotation * covarianceMatrix * qprotation.transpose();

  /// Make a unit vector for transformation matrix
  const double uPx = px / p;
  const double uPy = py / p;
  const double uPz = pz / p;
  const double uP = sqrt(uPx * uPx + uPy * uPy + uPz * uPz);

  const double cosphi = uPx / uP;
  const double sinphi = uPy / uP;
  const double rho = sqrt(uP * uP - uPz * uPz);
  const double sintheta = rho / uP;
  const double costheta = uPz / uP;

  Acts::BoundSymMatrix rotation = Acts::BoundSymMatrix::Zero();
  /// Rotate from d0, z0 to x,y,z. Again, assumes reference point of 0,0,0
  rotation(0,0) = cos(phiPos);
  rotation(1,0) = sin(phiPos);

  /// z pos trivial
  rotation(2,1) = 1;

  /// Leave uP in for clarity even though it is trivially 1
  rotation(3,2) = -1 * uP * sintheta * sinphi;
  rotation(3,3) = uP * costheta * cosphi;
  rotation(3,4) = sintheta * cosphi;
  rotation(4,2) = uP * sintheta * cosphi;
  rotation(4,3) = uP * costheta * sinphi;
  rotation(4,4) = sintheta * sinphi;
  rotation(5,3) = -1 * uP * sintheta;
  rotation(5,4) = costheta;
  
  printMatrix("Rotation matrix is ", rotation);

  matrix = rotation * covarianceMatrix * rotation.transpose();

  printMatrix("Rotated matrix is : ", matrix);

 

  return matrix;
}


void ActsCovarianceRotater::printMatrix(const std::string message,
					Acts::BoundSymMatrix matrix)
{
 
  if(m_verbosity)
    {
      std::cout << std::endl << message.c_str() << std::endl;
      for(int i = 0 ; i < matrix.rows(); ++i)
	{
	  std::cout << std::endl;
	  for(int j = 0; j < matrix.cols(); ++j)
	    {
	      std::cout << matrix(i,j) << ", ";
	    }
	}
    }

}
