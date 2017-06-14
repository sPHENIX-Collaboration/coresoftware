/*!
 *  \file		Measurement.h
 *  \brief		Measurement is the base class for input of the fitter.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef __PHGenFit_Measurement__
#define __PHGenFit_Measurement__

#include <climits>

#include <GenFit/AbsMeasurement.h>

namespace PHGenFit {

class Measurement {
public:
	//!ctor
	Measurement() :
			_measurement(NULL),
			_clusterID(UINT_MAX)
			{
	}
	;

	//!dtor
	~Measurement();

	//!
	genfit::AbsMeasurement* getMeasurement() {
		return _measurement;
	}

	unsigned int get_cluster_ID() const {
		return _clusterID;
	}

	void set_cluster_ID(unsigned int clusterId) {
		_clusterID = clusterId;
	}

protected:

	genfit::AbsMeasurement* _measurement;
	unsigned int _clusterID;

};
} //End of PHGenFit namespace

#endif //__PHGenFit_Measurement__
