/*!
 *  \file		SpacepointMeasurement.cc
 *  \brief		Handles the palnar type of measurements.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "SpacepointMeasurement.h"

#include <GenFit/SpacepointMeasurement.h>

#include <TMatrixDSymfwd.h>                // for TMatrixDSym
#include <TMatrixTSym.h>
#include <TVector3.h>
#include <TVectorDfwd.h>
#include <TVectorT.h>

namespace PHGenFit
{
void SpacepointMeasurement::init(const TVector3& pos, const TMatrixDSym& cov)
{
  int nDim = 3;
  TVectorD hitCoords(nDim);
  TMatrixDSym hitCov(nDim);

  hitCoords(0) = pos.X();
  hitCoords(1) = pos.Y();
  hitCoords(2) = pos.Z();

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      hitCov(i, j) = cov(i, j);

  int measurementCounter_ = 0;
  _measurement = new genfit::SpacepointMeasurement(hitCoords, hitCov, -1,
                                                   measurementCounter_,
                                                   nullptr);
}

SpacepointMeasurement::SpacepointMeasurement(const TVector3& pos, const double resolution)
{
  TMatrixDSym cov(3);
  cov.Zero();
  cov(0, 0) = resolution * resolution;
  cov(1, 1) = resolution * resolution;
  cov(2, 2) = resolution * resolution;
  init(pos, cov);
}

SpacepointMeasurement::SpacepointMeasurement(const TVector3& pos, const TVector3& resolution)
{
  TMatrixDSym cov(3);
  cov.Zero();
  cov(0, 0) = resolution.X() * resolution.X();
  cov(1, 1) = resolution.Y() * resolution.Y();
  cov(2, 2) = resolution.Z() * resolution.Z();
  init(pos, cov);
}
  
/*!
 * Ctor
 * \param pos measurement position
 * \param covariance matrix
 */
SpacepointMeasurement::SpacepointMeasurement(const TVector3& pos, const TMatrixDSym& cov)
{
  init(pos, cov);
}

}  // namespace PHGenFit
