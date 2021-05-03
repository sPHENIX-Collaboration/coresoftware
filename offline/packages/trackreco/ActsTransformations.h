#ifndef TRACKRECO_ACTSTRANSFORMATIONS_H
#define TRACKRECO_ACTSTRANSFORMATIONS_H

#include <trackbase/TrkrDefs.h>

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <trackbase_historic/SvtxTrack.h>

#include <ActsExamples/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>
#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

#include <boost/bimap.hpp>

/// std (and the like) includes
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

using SourceLink = ActsExamples::TrkrClusterSourceLink;

using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;
using Measurement = Acts::Measurement<ActsExamples::TrkrClusterSourceLink,
                                      Acts::BoundIndices,
                                      Acts::eBoundLoc0,
                                      Acts::eBoundLoc1>;

typedef boost::bimap<TrkrDefs::cluskey, unsigned int> CluskeyBimap;


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
 ActsTransformations()
   : m_verbosity(false)
  {}
  virtual ~ActsTransformations(){}
  
  /// Rotates an SvtxTrack covariance matrix from (x,y,z,px,py,pz) global
  /// cartesian coordinates to (d0, z0, phi, theta, q/p, time) coordinates for
  /// Acts. The track fitter performs the fitting with respect to the nominal
  /// origin of sPHENIX, so we rotate accordingly
  Acts::BoundSymMatrix rotateSvtxTrackCovToActs(const SvtxTrack *track,
						Acts::GeometryContext geoCtxt);
  
  /// Same as above, but rotate from Acts basis to global (x,y,z,px,py,pz)
  Acts::BoundSymMatrix rotateActsCovToSvtxTrack(
                       const Acts::BoundTrackParameters params,
		       Acts::GeometryContext geoCtxt);

  void setVerbosity(int verbosity) {m_verbosity = verbosity;}

  void printMatrix(const std::string &message, Acts::BoundSymMatrix matrix);

  /// Calculate the DCA for a given Acts fitted track parameters and 
  /// vertex
  void calculateDCA(const Acts::BoundTrackParameters param,
		    Acts::Vector3D vertex,
		    Acts::BoundSymMatrix cov,
		    Acts::GeometryContext geoCtxt,
		    float &dca3Dxy,
		    float &dca3Dz,
		    float &dca3DxyCov,
		    float &dca3DzCov);

  void fillSvtxTrackStates(const Trajectory traj, 
			   const size_t &trackTip,
			   SvtxTrack *svtxTrack,
			   Acts::GeometryContext geoContext);

 private:
  int m_verbosity;


};


#endif
