/*!
 *  \file		PlanarMeasurement.h
 *  \brief		Handles the palnar type of measurements.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef PHGENFIT_PLANARMEASUREMENT_H
#define PHGENFIT_PLANARMEASUREMENT_H

#include "Measurement.h"

class TVector3;

namespace PHGenFit
{
class PlanarMeasurement : public Measurement
{
 public:
  //!ctor
  PlanarMeasurement(const TVector3& pos, const TVector3& u, const TVector3& v, const double du, const double dv);

  PlanarMeasurement(const TVector3& pos, const TVector3& n, const double du, const double dv);

  void init(const TVector3& pos, const TVector3& u, const TVector3& v, const double du, const double dv);

  //!dtor
  ~PlanarMeasurement() override {}

 protected:
};
}  // namespace PHGenFit

#endif
