#ifndef TRACKRECO_ACTSALIGNMENTSTATES_H
#define TRACKRECO_ACTSALIGNMENTSTATES_H

#include <trackbase/TrkrDefs.h>

#include <Acts/Definitions/Algebra.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase_historic/SvtxAlignmentState.h>

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
  void clusterVersion(const int v) { m_clusterVersion = v; }
  void fillAlignmentStateMap(Trajectory traj,
                             SvtxTrack* track);
  void verbosity(const int verb) { m_verbosity = verb; }
  void analyticGlobalDer(bool a) { m_analytic = a; }
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
  
 private:
  void makeTpcGlobalCorrections(TrkrDefs::cluskey cluster_key, short int crossing, Acts::Vector3& global);
  std::vector<Acts::Vector3> getDerivativesAlignmentAngles(Acts::Vector3& global, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster, Surface surface, int crossing);
  void getGlobalDerivatives(std::vector<Acts::Vector3>& anglederivs, SvtxAlignmentState::GlobalMatrix& analytic);
  Acts::Transform3 makePerturbationTransformation(Acts::Vector3 angles);
  float convertTimeToZ(TrkrDefs::cluskey cluster_key, TrkrCluster* cluster);

  bool m_analytic = true;
  int m_verbosity = 0;
  int m_clusterVersion = 4;

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
};

#endif
