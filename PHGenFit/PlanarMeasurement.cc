/*!
 *  \file		PlanarMeasurement.cc
 *  \brief		Handles the palnar type of measurements.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include <GenFit/PlanarMeasurement.h>

#include "PlanarMeasurement.h"


namespace PHGenFit {


PlanarMeasurement::PlanarMeasurement(TVector3 pos, TVector3 u, TVector3 v, double du, double dv)
{

	int nDim = 2;
	TVectorD hitCoords(nDim);
	TMatrixDSym hitCov(nDim);

	hitCoords(0) = 0;
	hitCoords(1) = 0;

	hitCov(0,0) = du*du;
	hitCov(1,1) = dv*dv;

	genfit::SharedPlanePtr plane(
			new genfit::DetPlane(pos, u, v));

	int measurementCounter_ = 0;
	_measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, -1,
									measurementCounter_,
									nullptr);

	static_cast<genfit::PlanarMeasurement*>(_measurement)->setPlane(
			plane, measurementCounter_);
}

PlanarMeasurement::PlanarMeasurement(TVector3 pos, TVector3 n, double du, double dv)
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
	PlanarMeasurement(pos, u, v, du, dv);
}

PlanarMeasurement::~PlanarMeasurement()
{
	delete _measurement;
}



} //End of PHGenFit namespace
