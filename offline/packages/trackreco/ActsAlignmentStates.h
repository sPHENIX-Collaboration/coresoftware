#ifndef TRACKRECO_ACTSALIGNMENTSTATES_H
#define TRACKRECO_ACTSALIGNMENTSTATES_H

#include <trackbase/TrkrDefs.h>

#include <Acts/Definitions/Algebra.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase_historic/SvtxAlignmentState.h>

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

  ActsAlignmentStates() {}
  ~ActsAlignmentStates() {}
  
  void fillAlignmentStateMap(const ActsTrackFittingAlgorithm::TrackContainer& tracks,
			     const std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
                             SvtxTrack* track,
                             const ActsTrackFittingAlgorithm::MeasurementContainer& measurements);
  void verbosity(const int verb) { m_verbosity = verb; }
  void distortionContainers(TpcDistortionCorrectionContainer* stat,
                            TpcDistortionCorrectionContainer* average,
                            TpcDistortionCorrectionContainer* fluc)
  {
    m_dcc_static = stat;
    m_dcc_average = average;
    m_dcc_fluctuation = fluc;
  }
  void actsGeometry(ActsGeometry* geom) { m_tGeometry = geom; }
  void clusters(TrkrClusterContainer* clus) { m_clusterMap = clus; }
  void stateMap(SvtxAlignmentStateMap* map) { m_alignmentStateMap = map; }
  void fieldMap(std::string& fieldmap) {m_fieldMap = fieldmap; }

 private:
  void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global);

  std::pair<Acts::Vector3, Acts::Vector3> get_projectionXY(const Acts::Surface& surface, const Acts::Vector3& tangent);
  SvtxAlignmentState::GlobalMatrix makeGlobalDerivatives(const Acts::Vector3& OM, const std::pair<Acts::Vector3, Acts::Vector3>& projxy);

  int m_verbosity = 0;

  float sensorAngles[3] = {0.1, 0.1, 0.2};  // perturbation values for each alignment angle

  ClusterErrorPara m_clusErrPara;
  TpcDistortionCorrection m_distortionCorrection;
  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;

  SvtxAlignmentStateMap* m_alignmentStateMap = nullptr;
  TrkrClusterContainer* m_clusterMap = nullptr;
  ActsGeometry* m_tGeometry = nullptr;
  TpcDistortionCorrectionContainer* m_dcc_static = nullptr;
  TpcDistortionCorrectionContainer* m_dcc_average = nullptr;
  TpcDistortionCorrectionContainer* m_dcc_fluctuation = nullptr;
  std::string m_fieldMap = "";
};

#endif
