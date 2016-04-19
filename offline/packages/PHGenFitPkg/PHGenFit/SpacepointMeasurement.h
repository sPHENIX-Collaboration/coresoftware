/*!
 *  \file		SpacepointMeasurement.h
 *  \brief		Handles the palnar type of measurements.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */


#ifndef __PHGenFit_SpacepointMeasurement__
#define __PHGenFit_SpacepointMeasurement__

#include "Measurement.h"

class TVector3;

namespace PHGenFit {

class SpacepointMeasurement : public Measurement
{
public:
	//!ctor
	SpacepointMeasurement(const TVector3& pos, const double resolution);

	void init(const TVector3& pos, const double resolution);

	//!dtor
	~SpacepointMeasurement();

protected:

	};
} //End of PHGenFit namespace

#endif //__PHGenFit_SpacepointMeasurement__
