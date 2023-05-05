// $Id: $

/*!
 * \file RawClusterUtility.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef CALOBASE_RAWCLUSTERUTILITY_H
#define CALOBASE_RAWCLUSTERUTILITY_H

#include "RawCluster.h"

#include <CLHEP/Vector/ThreeVector.h>

/*!
 * \brief RawClusterUtility
 */
class RawClusterUtility
{
 public:
  virtual ~RawClusterUtility() {}

  //
  //! polar angle of cluster with respect to given vertex position
  static inline float
  GetPolarAngle(const RawCluster& cluster, const CLHEP::Hep3Vector vertex) { return (cluster.get_position() - vertex).getTheta(); }
  //
  //! azimuth angle of cluster with respect to given vertex position
  static inline float
  GetAzimuthAngle(const RawCluster& cluster, const CLHEP::Hep3Vector vertex) { return (cluster.get_position() - vertex).getPhi(); }
  //
  //! Pseudorapidity of cluster with respect to given vertex position
  static inline float
  GetPseudorapidity(const RawCluster& cluster, const CLHEP::Hep3Vector vertex) { return (cluster.get_position() - vertex).pseudoRapidity(); }
  //
  //! Transverse energy of cluster with respect to given vertex position
  static inline float
  GetET(const RawCluster& cluster, const CLHEP::Hep3Vector vertex)
  {
    const CLHEP::Hep3Vector cluster_vec(cluster.get_position() - vertex);
    const double mag = cluster_vec.mag();

    return mag > 0 ? cluster_vec.perp() / cluster_vec.mag() * cluster.get_energy() : 0;
  }
  //
  //! Transverse core energy of cluster with respect to given vertex position
  static inline float
  GetETCore(const RawCluster& cluster, const CLHEP::Hep3Vector vertex)
  {
    const CLHEP::Hep3Vector cluster_vec(cluster.get_position() - vertex);

    return cluster_vec.perp() / cluster_vec.mag() * cluster.get_ecore();
  }
  //
  //! Get energy vector of cluster based on cluster energy
  static inline CLHEP::Hep3Vector GetEVec(const RawCluster& cluster, const CLHEP::Hep3Vector vertex) { return (cluster.get_position() - vertex).unit() * cluster.get_energy(); }
  //! Get energy vector of cluster based on E_core
  static inline CLHEP::Hep3Vector GetECoreVec(const RawCluster& cluster, const CLHEP::Hep3Vector vertex) { return (cluster.get_position() - vertex).unit() * cluster.get_ecore(); }

 private:
  //! not intended to make an instance
  RawClusterUtility() {}
};

#endif
