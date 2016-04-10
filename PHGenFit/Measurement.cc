/*!
 *  \file		Measurement.cc
 *  \brief		Measurement is the base class for input of the fitter.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#include "Measurement.h"

namespace PHGenFit {

Measurement::~Measurement()
{
	delete _measurement;
}
} //End of PHGenFit namespace
