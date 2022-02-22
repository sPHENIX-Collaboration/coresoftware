/*!
 *  \file		Measurement.h
 *  \brief		Measurement is the base class for input of the fitter.
 *  \author		Haiwang Yu <yuhw@nmsu.edu>
 */

#ifndef PHGENFIT_MEASUREMENT_H
#define PHGENFIT_MEASUREMENT_H

#include <GenFit/AbsMeasurement.h>
#include <trackbase/TrkrDefs.h>

#include <climits>

namespace PHGenFit
{
class Measurement
{
 public:
  //!ctor
  Measurement()
    : _measurement(NULL)
    , _clusterID(UINT_MAX){};

  //!dtor
  virtual ~Measurement() {}

  //!
  genfit::AbsMeasurement* getMeasurement()
  {
    return _measurement;
  }

  // old tracking
  unsigned int get_cluster_ID() const
  {
    return _clusterID;
  }
  void set_cluster_ID(unsigned int clusterId)
  {
    _clusterID = clusterId;
  }

  // new tracking
  TrkrDefs::cluskey get_cluster_key() const
  {
    return _clusterkey;
  }
  void set_cluster_key(TrkrDefs::cluskey clusterkey)
  {
    _clusterkey = clusterkey;
  }

 protected:
  genfit::AbsMeasurement* _measurement;
  unsigned int _clusterID;
  TrkrDefs::cluskey _clusterkey;

};
}  // namespace PHGenFit

#endif
