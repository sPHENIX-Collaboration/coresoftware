// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHACTSGSF_H
#define PHACTSGSF_H

#include <fun4all/SubsysReco.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/ActsSourceLink.h>
#include <trackbase/Calibrator.h>
#include <trackbase/ClusterErrorPara.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

#include <trackbase/ActsTrackFittingAlgorithm.h>

#include <string>

class PHCompositeNode;
class ActsGeometry;
class TrkrClusterContainer;
class SvtxTrackMap;
class SvtxVertexMap;
class SvtxTrack;

using SourceLink = ActsSourceLink;
using FitResult = ActsTrackFittingAlgorithm::TrackFitterResult;
using Trajectory = ActsExamples::Trajectories;
using Measurement = Acts::Measurement<Acts::BoundIndices, 2>;
using SurfacePtrVec = std::vector<const Acts::Surface*>;
using SourceLinkVec = std::vector<Acts::SourceLink>;

class PHActsGSF : public SubsysReco
{
 public:
  PHActsGSF(const std::string& name = "PHActsGSF");

  ~PHActsGSF() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_pp_mode(bool mode) {m_pp_mode = mode;}

 private:
  int getNodes(PHCompositeNode* topNode);
  std::shared_ptr<Acts::PerigeeSurface> makePerigee(SvtxTrack* track) const;
  ActsTrackFittingAlgorithm::TrackParameters makeSeed(
      SvtxTrack* track,
      std::shared_ptr<Acts::PerigeeSurface> psurf) const;
  //  SourceLinkVec getSourceLinks(TrackSeed* track,
  //                         ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
  //                         const short int& crossing);
  ActsTrackFittingAlgorithm::TrackFitterResult fitTrack(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const ActsTrackFittingAlgorithm::TrackParameters& seed,
      const ActsTrackFittingAlgorithm::GeneralFitterOptions& options,
      const CalibratorAdapter& calibrator,
      ActsTrackFittingAlgorithm::TrackContainer& tracks);

  void updateTrack(FitResult& result, SvtxTrack* track,
                   ActsTrackFittingAlgorithm::TrackContainer& tracks);
  void updateSvtxTrack(std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
                       Trajectory::IndexedParameters& paramsMap,
                       ActsTrackFittingAlgorithm::TrackContainer& tracks,
                       SvtxTrack* track);
  ActsGeometry* m_tGeometry = nullptr;
  TrkrClusterContainer* m_clusterContainer = nullptr;
  SvtxTrackMap* m_trackMap = nullptr;
  SvtxVertexMap* m_vertexMap = nullptr;
  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;

  alignmentTransformationContainer* m_alignmentTransformationMap = nullptr;  // added for testing purposes
  alignmentTransformationContainer* m_alignmentTransformationMapTransient = nullptr;  
  std::set< Acts::GeometryIdentifier> m_transient_id_set;
  Acts::GeometryContext m_transient_geocontext;

  TpcDistortionCorrectionContainer* m_dccStatic = nullptr;
  TpcDistortionCorrectionContainer* m_dccAverage = nullptr;
  TpcDistortionCorrectionContainer* m_dccFluctuation{nullptr};
  TpcDistortionCorrection m_distortionCorrection;

  std::string m_trackMapName = "SvtxTrackMap";
  unsigned int m_pHypothesis = 11;

  bool m_pp_mode = false;

  ClusterErrorPara _ClusErrPara;

  ActsTrackFittingAlgorithm::Config m_fitCfg;
};

#endif  // PHACTSGSF_H
