/*!
 *  \file		PlanarMeasurement.cc
 *  \brief		Handles the palnar type of measurements.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "PlanarMeasurement.h"

#include <GenFit/PlanarMeasurement.h>
#include <GenFit/DetPlane.h>
#include <GenFit/SharedPlanePtr.h>     // for SharedPlanePtr
#include <GenFit/StateOnPlane.h>

#include <TMatrixDSymfwd.h>            // for TMatrixDSym
#include <TMatrixTSym.h>               // for TMatrixTSym
#include <TVectorDfwd.h>               // for TVectorD
#include <TVectorT.h>                  // for TVectorT
#include <TVector3.h>

namespace PHGenFit
{
void PlanarMeasurement::init(const TVector3& pos, const TVector3& u, const TVector3& v, const double du, const double dv)
{
  int nDim = 2;
  TVectorD hitCoords(nDim);
  TMatrixDSym hitCov(nDim);

  hitCoords(0) = 0;
  hitCoords(1) = 0;

  hitCov(0, 0) = du * du;
  hitCov(1, 1) = dv * dv;

  genfit::SharedPlanePtr plane(
      new genfit::DetPlane(pos, u, v));

  int measurementCounter_ = 0;
  _measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, -1,
                                               measurementCounter_,
                                               nullptr);

  static_cast<genfit::PlanarMeasurement*>(_measurement)->setPlane(plane, measurementCounter_);
}

PlanarMeasurement::PlanarMeasurement(const TVector3& pos, const TVector3& u, const TVector3& v, const double du, const double dv)
{
  init(pos, u, v, du, dv);
}

PlanarMeasurement::PlanarMeasurement(const TVector3& pos, const TVector3& n, const double du, const double dv)
{
  /*!
	 *  Mainly for "n" in xy plane, or n.z() = 0;
	 *  In this case, According to https://root.cern.ch/doc/master/TVector3_8h_source.html#l00301
	 *  u = n.Orthogonal() is always in xy plane
	 *  v = n x u, is along z direction.
	 *  If z is not 0, but z <= x or y, then u will still be in xy plane, phi direction error, but v will not be along z axis.
	 */
  TVector3 u = n.Orthogonal();
  TVector3 v = n.Cross(u);
  init(pos, u, v, du, dv);
}

}  // namespace PHGenFit
