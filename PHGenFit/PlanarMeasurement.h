/*!
 *  \file		PlanarMeasurement.h
 *  \brief		Handles the palnar type of measurements.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */


#ifndef __PHGenFit_PlanarPlanarMeasurement__
#define __PHGenFit_PlanarPlanarMeasurement__

#include "Measurement.h"

class TVector3;

namespace PHGenFit {

class PlanarMeasurement : public Measurement
{
public:
	//!ctor
	PlanarMeasurement(TVector3 pos, TVector3 u, TVector3 v, double du, double dv);

	//!dtor
	~PlanarMeasurement();

protected:

	};
} //End of PHGenFit namespace

#endif //__PHGenFit_PlanarPlanarMeasurement__
