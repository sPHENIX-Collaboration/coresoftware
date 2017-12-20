// $Id: $

/*!
 * \file RawClusterUtility.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4CEMC_RAWCLUSTERUTILITY_H_
#define SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4CEMC_RAWCLUSTERUTILITY_H_

#include <CLHEP/Vector/ThreeVector.h>

class RawCluster;

/*!
 * \brief RawClusterUtility
 */
class RawClusterUtility
{
  virtual ~RawClusterUtility();

  //
  //! polar angle of cluster with respect to given vertex position
  static float
  GetPolarAngle(const RawCluster& cluster, const CLHEP::Hep3Vector vertex) { return (cluster.get_position() - vertex).getTheta(); }
  //
  //! azimuth angle of cluster with respect to given vertex position
  static float
  GetAzimuthAngle(const RawCluster& cluster, const CLHEP::Hep3Vector vertex) { return (cluster.get_position() - vertex).getPhi(); }
  //
  //! Pseudorapidity of cluster with respect to given vertex position
  static float
  GetPseudorapidity(const RawCluster& cluster, const CLHEP::Hep3Vector vertex) { return (cluster.get_position() - vertex).pseudoRapidity(); }
  //
  //! Transverse energy of cluster with respect to given vertex position
  static float
  GetET(const RawCluster& cluster, const CLHEP::Hep3Vector vertex)
  {
    const CLHEP::Hep3Vector cluster_vec(cluster.get_position() - vertex);

    return cluster_vec.perp() / cluster_vec.mag() * cluster.get_energy();
  }
  //
  //! Transverse core energy of cluster with respect to given vertex position
  static float
  GetETCore(const RawCluster& cluster, const CLHEP::Hep3Vector vertex)
  {
    const CLHEP::Hep3Vector cluster_vec(cluster.get_position() - vertex);

    return cluster_vec.perp() / cluster_vec.mag() * cluster.get_ecore();
  }

 private:
  //! not intended to make an instance
  RawClusterUtility();
};

#endif /* SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_G4CEMC_RAWCLUSTERUTILITY_H_ */
