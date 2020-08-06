#ifndef TRACKRECO_ACTSCOVARIANCEROTATER_H
#define TRACKRECO_ACTSCOVARIANCEROTATER_H

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <trackbase_historic/SvtxTrack.h>

#include <ACTFW/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ACTFW/EventData/TrkrClusterSourceLink.hpp>

/// std (and the like) includes
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

using SourceLink = FW::Data::TrkrClusterSourceLink;


/**
 * This is a helper class for rotating track covariance matrices to and from
 * the basis that Acts expects. The covariance matrix is nominally given in the
 * global basis (x,y,z,px,py,pz). Acts expects the covariance matrix in a local
 * basis with respect to the given reference point that is provided as an
 * option to the KalmanFitter. 
 */
class ActsCovarianceRotater
{
  public:
  ActsCovarianceRotater()
    : m_verbosity(false)
    {}
  virtual ~ActsCovarianceRotater(){}
  
  /// Rotates an SvtxTrack covariance matrix from (x,y,z,px,py,pz) global
  /// cartesian coordinates to (d0, z0, phi, theta, q/p, time) coordinates for
  /// Acts. The track fitter performs the fitting with respect to the nominal
  /// origin of sPHENIX, so we rotate accordingly
  Acts::BoundSymMatrix rotateSvtxTrackCovToActs(const SvtxTrack *track);
  
  /// Same as above, but rotate from Acts basis to global (x,y,z,px,py,pz)
  Acts::BoundSymMatrix rotateActsCovToSvtxTrack(const Acts::BoundParameters params);

  void setVerbosity(int verbosity) {m_verbosity = verbosity;}

  void printMatrix(const std::string &message, Acts::BoundSymMatrix matrix);

 private:
  int m_verbosity;


};


#endif
