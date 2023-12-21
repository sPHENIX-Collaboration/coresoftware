#ifndef TRACKRECO_MAKESOURCELINKS_H
#define TRACKRECO_MAKESOURCELINKS_H

#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase/alignmentTransformationContainer.h>
#include <trackbase/ClusterErrorPara.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Definitions/Algebra.hpp>

#include <trackbase_historic/SvtxTrack.h>

/// std (and the like) includes
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

using SourceLink = ActsSourceLink;
using SourceLinkVec = std::vector<Acts::SourceLink>;

// forward declarations
class SvtxTrack;
class SvtxTrackState;
class TrkrCluster;
class TrkrClusterContainer;
class TpcDistortionCorrectionContainer;
class TrackSeed;

/**
   This class contains the code that generates Acts source links and makes the modified GeoContext 
   that applies cluster crossing and distortion corrections in the alignment transformations.
   It is used by all of the Acts track fitting modules.
 */
class MakeSourceLinks
{
 public:
 MakeSourceLinks() = default;
    
   void setVerbosity(int verbosity) {m_verbosity = verbosity;}

 void set_pp_mode(bool ispp) { m_pp_mode = ispp; }

  void ignoreLayer(int layer) { m_ignoreLayer.insert(layer); }
  
  SourceLinkVec getSourceLinks(TrackSeed* track,
			       ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
			       TrkrClusterContainer*  cluster_container,
			       ActsGeometry* tGeometry,
			       alignmentTransformationContainer* transformMapTransient,
			       std::set< Acts::GeometryIdentifier> transient_id_set,
			       short int crossing);

  void resetTransientTransformMap(
						   alignmentTransformationContainer* transformMapTransient,
						   std::set< Acts::GeometryIdentifier>& transient_id_set,
						   ActsGeometry* tGeometry );

 private:
  int m_verbosity = 0;
  bool m_pp_mode = false;  
  std::set<int> m_ignoreLayer;

  TpcClusterZCrossingCorrection _clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};
  
  /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;
  
  ClusterErrorPara _ClusErrPara;
  

};


#endif
