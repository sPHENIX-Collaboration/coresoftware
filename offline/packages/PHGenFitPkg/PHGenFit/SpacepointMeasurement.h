/*!
 *  \file		SpacepointMeasurement.h
 *  \brief		Handles the palnar type of measurements.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef PHGENFIT_SPACEPOINTMEASUREMENT_H
#define PHGENFIT_SPACEPOINTMEASUREMENT_H

#include "Measurement.h"

#include <TMatrixDSymfwd.h>

class TVector3;

namespace PHGenFit
{
class SpacepointMeasurement : public Measurement
{
 public:
  /*!
	 * Ctor
	 * \param pos measurement position
	 * \param resolution standard dev for diagnal elements of the cov, other elements are zero
	 */
  SpacepointMeasurement(const TVector3& pos, const double resolution);
  /*!
	 * Ctor
	 * \param pos measurement position
	 * \param resolution standard dev for each diagnal element of the cov, other elements are zero
	 */
  SpacepointMeasurement(const TVector3& pos, const TVector3& resolution);
  /*!
	 * Ctor
	 * \param pos measurement position
	 * \param covariance matrix
	 */
  SpacepointMeasurement(const TVector3& pos, const TMatrixDSym& cov);

  void init(const TVector3& pos, const TMatrixDSym& cov);

  //!dtor
  ~SpacepointMeasurement() override {}

 protected:
};
}  // namespace PHGenFit

#endif
