#ifndef TRACKRECO_ACTSTRANSFORMATIONS_H
#define TRACKRECO_ACTSTRANSFORMATIONS_H

#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include "SvtxTrack.h"

#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

/// std (and the like) includes
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

class TrkrCluster;

using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;


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
  Acts::BoundSymMatrix rotateSvtxTrackCovToActs(const SvtxTrack *track,
						Acts::GeometryContext geoCtxt) const;
  
  /// Same as above, but rotate from Acts basis to global (x,y,z,px,py,pz)
  Acts::BoundSymMatrix rotateActsCovToSvtxTrack(
                       const Acts::BoundTrackParameters params,
		       Acts::GeometryContext geoCtxt) const;

  void setVerbosity(int verbosity) {m_verbosity = verbosity;}

  void printMatrix(const std::string &message, Acts::BoundSymMatrix matrix) const;

  /// Calculate the DCA for a given Acts fitted track parameters and 
  /// vertex
  void calculateDCA(const Acts::BoundTrackParameters param,
		    Acts::Vector3D vertex,
		    Acts::BoundSymMatrix cov,
		    Acts::GeometryContext geoCtxt,
		    float &dca3Dxy,
		    float &dca3Dz,
		    float &dca3DxyCov,
		    float &dca3DzCov) const;

  void fillSvtxTrackStates(const Trajectory& traj, 
			   const size_t &trackTip,
			   SvtxTrack *svtxTrack,
			   Acts::GeometryContext geoContext) const;
  
  Acts::Vector3F getGlobalPositionF(TrkrCluster* cluster,
				    ActsSurfaceMaps* surfMaps,
				    ActsTrackingGeometry *tGeometry) const;
  Acts::Vector3D getGlobalPosition(TrkrCluster* cluster,
				   ActsSurfaceMaps* surfMaps,
				   ActsTrackingGeometry *tGeometry) const;
  Surface getSurface(TrkrCluster* cluster,
		     ActsSurfaceMaps* surfMaps) const;
  
 private:
  int m_verbosity = 0;

  Surface getSiliconSurface(TrkrDefs::hitsetkey hitsetkey,
			    ActsSurfaceMaps *maps) const;
  Surface getTpcSurface(TrkrDefs::hitsetkey hitsetkey,
			TrkrDefs::subsurfkey surfkey,
			ActsSurfaceMaps *maps) const;
  Surface getMMSurface(TrkrDefs::hitsetkey hitsetkey,
		       ActsSurfaceMaps *maps) const;
  

};


#endif
