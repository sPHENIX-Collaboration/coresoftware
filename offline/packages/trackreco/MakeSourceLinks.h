#ifndef TRACKRECO_MAKESOURCELINKS_H
#define TRACKRECO_MAKESOURCELINKS_H

#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase/alignmentTransformationContainer.h>
#include <trackbase/ClusterErrorPara.h>

#include <tpc/TpcClusterMover.h>

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
class TpcGlobalPositionWrapper;
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

 void initialize(PHG4TpcCylinderGeomContainer* cellgeo);

  void setVerbosity(int verbosity) {m_verbosity = verbosity;}

 void set_pp_mode(bool ispp) { m_pp_mode = ispp; }

  void ignoreLayer(int layer) { m_ignoreLayer.insert(layer); }

  SourceLinkVec getSourceLinks(
    TrackSeed* /*seed*/,
    ActsTrackFittingAlgorithm::MeasurementContainer& /*measurements*/,
    TrkrClusterContainer* /*clusters*/,
    ActsGeometry* /*geometry*/,
    const TpcGlobalPositionWrapper& /*globalpositionWrapper*/,
    alignmentTransformationContainer* /*transformMapTransient*/,
    std::set< Acts::GeometryIdentifier>& /*transient_id_set*/,
    short int /*crossing*/);

  void resetTransientTransformMap(
    alignmentTransformationContainer* /*transformMapTransient*/,
    std::set< Acts::GeometryIdentifier>& /*transient_id_set*/,
    ActsGeometry* /*tGeometry*/ );

  SourceLinkVec getSourceLinksClusterMover(
    TrackSeed* /*seed*/,
    ActsTrackFittingAlgorithm::MeasurementContainer& /*measurements*/,
    TrkrClusterContainer* /*clusters*/,
    ActsGeometry* /*geometry*/,
    const TpcGlobalPositionWrapper& /*globalpositionWrapper*/,
    short int crossing
    );

 private:
  int m_verbosity = 0;
  bool m_pp_mode = false;
  std::set<int> m_ignoreLayer;

  TpcClusterMover _clusterMover;

  ClusterErrorPara _ClusErrPara;


};


#endif
