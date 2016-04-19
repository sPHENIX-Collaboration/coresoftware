/*!
 *  \file		SpacepointMeasurement.cc
 *  \brief		Handles the palnar type of measurements.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include <GenFit/SpacepointMeasurement.h>

#include "SpacepointMeasurement.h"


namespace PHGenFit {

void SpacepointMeasurement::init(const TVector3& pos, const double resolution)
{

	int nDim = 3;
	TVectorD hitCoords(nDim);
	TMatrixDSym hitCov(nDim);

	hitCoords(0) = pos.X();
	hitCoords(1) = pos.Y();
	hitCoords(2) = pos.Z();

	hitCov(0,0) = resolution*resolution;
	hitCov(1,1) = resolution*resolution;
	hitCov(2,2) = resolution*resolution;

	int measurementCounter_ = 0;
	_measurement = new genfit::SpacepointMeasurement(hitCoords, hitCov, -1,
									measurementCounter_,
									nullptr);

}

SpacepointMeasurement::SpacepointMeasurement(const TVector3& pos, const double resolution)
{
	init(pos, resolution);
}

SpacepointMeasurement::~SpacepointMeasurement()
{
	//delete _measurement;
}



} //End of PHGenFit namespace
