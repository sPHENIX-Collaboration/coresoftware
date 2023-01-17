#ifndef TRACKRECO_ACTSTRANSFORMATIONS_H
#define TRACKRECO_ACTSTRANSFORMATIONS_H

#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>

#include "SvtxTrack.h"

/// std (and the like) includes
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

// forward declarations
class SvtxTrack;
class SvtxTrackState;
class TrkrCluster;

/**
 * This is a helper class for rotating track covariance matrices to and from
 * the basis that Acts expects. The covariance matrix is nominally given in the
 * global basis (x,y,z,px,py,pz). Acts expects the covariance matrix in a local
 * basis with respect to the given reference point that is provided as an
 * option to the KalmanFitter. 
 */
class ActsTransformations
{
 public:
 ActsTransformations() = default;
    
  /// Rotates an SvtxTrack covariance matrix from (x,y,z,px,py,pz) global
  /// cartesian coordinates to (d0, z0, phi, theta, q/p, time) coordinates for
  /// Acts. The track fitter performs the fitting with respect to the nominal
  /// origin of sPHENIX, so we rotate accordingly
  Acts::BoundSymMatrix rotateSvtxTrackCovToActs(const SvtxTrack* ) const;
  
  /// Rotates an SvtxTrack state covariance matrix from (x,y,z,px,py,pz) global
  /// cartesian coordinates to (d0, z0, phi, theta, q/p, time) coordinates for
  /// Acts. The track fitter performs the fitting with respect to the nominal
  /// origin of sPHENIX, so we rotate accordingly
  Acts::BoundSymMatrix rotateSvtxTrackCovToActs(const SvtxTrackState* ) const;

  /// Rotates an Acts covariance matrix from (d0, z0, phi, theta, q/p, time) local curvilinear coordinates
  /// to global cartesian coordinates (x,y,z,px,py,pz) coordinates
  Acts::BoundSymMatrix rotateActsCovToSvtxTrack( const ActsTrackFittingAlgorithm::TrackParameters& ) const;

  void setVerbosity(int verbosity) {m_verbosity = verbosity;}

  void printMatrix(const std::string &message, const Acts::BoundSymMatrix& matrix) const;

  /// Calculate the DCA for a given Acts fitted track parameters and 
  /// vertex
  void calculateDCA(const ActsTrackFittingAlgorithm::TrackParameters param,
		    Acts::Vector3 vertex,
		    Acts::BoundSymMatrix cov,
		    Acts::GeometryContext& geoCtxt,
		    float &dca3Dxy,
		    float &dca3Dz,
		    float &dca3DxyCov,
		    float &dca3DzCov) const;

  void fillSvtxTrackStates(const Acts::MultiTrajectory<Acts::VectorMultiTrajectory>& traj, 
			   const size_t& trackTip,
			   SvtxTrack *svtxTrack,
			   Acts::GeometryContext& geoContext) const;
  
 private:
  int m_verbosity = 0;
  

};


#endif
