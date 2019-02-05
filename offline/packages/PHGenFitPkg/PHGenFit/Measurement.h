/*!
 *  \file		Measurement.h
 *  \brief		Measurement is the base class for input of the fitter.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef PHGENFIT_MEASUREMENT_H
#define PHGENFIT_MEASUREMENT_H


#include <GenFit/AbsMeasurement.h>

#include <climits>

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
			~Measurement(){}

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

#endif
