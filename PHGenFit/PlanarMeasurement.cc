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

PlanarMeasurement::~PlanarMeasurement()
{
	delete _measurement;
}



} //End of PHGenFit namespace
