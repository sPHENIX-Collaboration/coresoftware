// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHACTSGSF_H
#define PHACTSGSF_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrectionContainer.h>
#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcClusterMover.h>

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp>
#include <ActsExamples/EventData/Trajectories.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/EventData/IndexSourceLink.hpp>

class PHCompositeNode;
class ActsGeometry;
class TrkrClusterContainer;
class SvtxTrackMap;
class SvtxVertexMap;

using SourceLink = ActsExamples::IndexSourceLink;
using FitResult = Acts::KalmanFitterResult;
using Trajectory = ActsExamples::Trajectories;
using Measurement = Acts::Measurement<Acts::BoundIndices,2>;
using SurfacePtrVec = std::vector<const Acts::Surface*>;
using SourceLinkVec = std::vector<SourceLink>;

class PHActsGSF : public SubsysReco
{
 public:

  PHActsGSF(const std::string &name = "PHActsGSF");

  ~PHActsGSF() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:

  int getNodes(PHCompositeNode* topNode);
  std::shared_ptr<Acts::PerigeeSurface> makePerigee(SvtxTrack* track) const;
  ActsExamples::TrackParameters makeSeed(
		        SvtxTrack* track,
			std::shared_ptr<Acts::PerigeeSurface> psurf) const;
  SourceLinkVec getSourceLinks(TrackSeed* track, 
			       ActsExamples::MeasurementContainer& measurements,
			       const short int& crossing);
  ActsExamples::TrackFittingAlgorithm::TrackFitterResult fitTrack(
      const std::vector<std::reference_wrapper<const SourceLink>>& sourceLinks,
      const ActsExamples::TrackParameters& seed,
      const ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions& options);

  void updateTrack(const FitResult& result, SvtxTrack *track);
  void updateSvtxTrack(const Trajectory& traj, SvtxTrack* track);
  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;

  TpcDistortionCorrectionContainer* m_dccStatic = nullptr;
  TpcDistortionCorrectionContainer* m_dccAverage = nullptr;
  TpcDistortionCorrectionContainer* m_dccFluctuation{nullptr};
  TpcDistortionCorrection m_distortionCorrection;
  TpcClusterMover m_clusterMover;

  std::string m_trackMapName = "SvtxTrackMap";
  unsigned int m_pHypothesis = 11;
  
  ActsExamples::TrackFittingAlgorithm::Config m_fitCfg;
};

#endif // PHACTSGSF_H
