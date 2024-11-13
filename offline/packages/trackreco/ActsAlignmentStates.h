#ifndef TRACKRECO_ACTSALIGNMENTSTATES_H
#define TRACKRECO_ACTSALIGNMENTSTATES_H

#include <trackbase/TrkrDefs.h>

#include <Acts/Definitions/Algebra.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

#include <tpc/TpcClusterZCrossingCorrection.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase_historic/SvtxAlignmentState.h>

#include <tpc/TpcGlobalPositionWrapper.h>

#include <string>

class PHCompositeNode;
class SvtxTrack;
class SvtxAlignmentStateMap;
class ActsGeometry;
class TrkrClusterContainer;
/*
 * Helper class that contains functionality to fill alignment state
 * map from Acts track fitter
 */
class ActsAlignmentStates
{
 public:
  using Trajectory = ActsExamples::Trajectories;

  explicit ActsAlignmentStates() = default;

  void fillAlignmentStateMap(const ActsTrackFittingAlgorithm::TrackContainer& tracks,
    const std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
    SvtxTrack* track,
    const ActsTrackFittingAlgorithm::MeasurementContainer& measurements);

  //! load all relevant nodes
  void loadNodes( PHCompositeNode* /*topNode*/);

  //! verbosity
  void verbosity(const int verb)
  { m_verbosity = verb; }

  //! field map string
  void fieldMap(std::string& fieldmap)
  {m_fieldMap = fieldmap; }

 private:

  std::pair<Acts::Vector3, Acts::Vector3> get_projectionXY(const Acts::Surface& surface, const Acts::Vector3& tangent);

  SvtxAlignmentState::GlobalMatrix makeGlobalDerivatives(const Acts::Vector3& OM, const std::pair<Acts::Vector3, Acts::Vector3>& projxy);

  //! verbosity
  int m_verbosity = 0;

  //! global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  //! alignment state map
  SvtxAlignmentStateMap* m_alignmentStateMap = nullptr;

  //! cluster map
  TrkrClusterContainer* m_clusterMap = nullptr;

  //! acts geometry
  ActsGeometry* m_tGeometry = nullptr;

  //! field map
  std::string m_fieldMap;
};

#endif
