#ifndef TRACKRECO_PHACTSTRACKPROJECTION_H
#define TRACKRECO_PHACTSTRACKPROJECTION_H

#include "ActsPropagator.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>

#include <fun4all/SubsysReco.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Utilities/Result.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

#include <map>
#include <memory>
#include <string>

class PHCompositeNode;
class RawClusterContainer;
class RawTowerGeomContainer;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class TowerInfoContainer;
class TrkrClusterContainer;

/**
 * This class takes final fitted tracks from the Acts track fitting
 * and projects them out to cylinders with radius at the same radius
 * as the three calorimeters. Cluster matching is performed with the
 * projections and the SvtxTrack object is updated.
 */

class PHActsTrackProjection : public SubsysReco
{
 public:
  using BoundTrackParam =
      const Acts::BoundTrackParameters;
  using SurfacePtr = std::shared_ptr<const Acts::Surface>;
  using Trajectory = ActsExamples::Trajectories;
  using BoundTrackParamResult = ActsPropagator::BTPPairResult;

  PHActsTrackProjection(const std::string &name = "PHActsTrackProjection");

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void setConstFieldVal(float b) { m_constFieldVal = b; }

  /// Set an arbitrary radius to project to, in cm
  void setLayerRadius(SvtxTrack::CAL_LAYER layer, float rad)
  { m_caloRadii[layer] = rad; }

  // extrapolation mode
  enum class ExtrapolationMode
  {
    Legacy, // the default extrapolation mode, using fitter track parameters at origin
    Forward // uses the track state vector closest to the requested layer, before
  };

  /// extrapolation mode
  void setExtrapolationMode( const ExtrapolationMode value )
  { m_extrapolation_mode = value; }

 private:
  int getNodes(PHCompositeNode *topNode);
  int projectTracks(SvtxTrack::CAL_LAYER);

  // project a give track to all available calorimeter layers
  void projectTrack( SvtxTrack* ) const;

  /// Make Acts::CylinderSurface objects corresponding to the calos
  int makeCaloSurfacePtrs(PHCompositeNode *topNode);

  /// Update the SvtxTrack object with the track-cluster match
  void updateSvtxTrack(const ActsPropagator::BoundTrackParamPair &params,
                       SvtxTrack *svtxTrack,
                       SvtxTrack::CAL_LAYER) const;

  /// Objects containing the Acts track fit results
  ActsGeometry *m_tGeometry = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;

  /// cluster container
  TrkrClusterContainer *m_clusterContainer = nullptr;

  /// Objects to hold calorimeter information.
  std::map<SvtxTrack::CAL_LAYER, SurfacePtr> m_caloSurfaces;

  /// An optional map that allows projection to an arbitrary radius
  /// Results are written to the SvtxTrack based on the provided CAL_LAYER
  std::map<SvtxTrack::CAL_LAYER, float> m_caloRadii;

  /// constant field value
  float m_constFieldVal = 1.4;

  /// extrapolation mode
  ExtrapolationMode m_extrapolation_mode = ExtrapolationMode::Legacy;

};

#endif
