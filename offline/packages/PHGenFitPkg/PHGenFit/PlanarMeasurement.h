/*!
 *  \file		PlanarMeasurement.h
 *  \brief		Handles the palnar type of measurements.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */


#ifndef __PHGenFit_PlanarMeasurement__
#define __PHGenFit_PlanarMeasurement__

#include "Measurement.h"

class TVector3;

namespace PHGenFit {

class PlanarMeasurement : public Measurement
{
public:
	//!ctor
	PlanarMeasurement(const TVector3& pos, const TVector3& u, const TVector3& v, const double du, const double dv);

	PlanarMeasurement(const TVector3& pos, const TVector3& n, const double du, const double dv);

	void init(const TVector3& pos, const TVector3& u, const TVector3& v, const double du, const double dv);

	//!dtor
	~PlanarMeasurement();

protected:

	};
} //End of PHGenFit namespace

#endif //__PHGenFit_PlanarMeasurement__
