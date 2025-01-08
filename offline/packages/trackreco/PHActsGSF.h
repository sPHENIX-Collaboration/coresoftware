// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHACTSGSF_H
#define PHACTSGSF_H

#include "ActsEvaluator.h"

#include <fun4all/SubsysReco.h>

#include <tpc/TpcGlobalPositionWrapper.h>

#include <trackbase/ActsSourceLink.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase/Calibrator.h>
#include <trackbase/ClusterErrorPara.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>


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

  //! constructor
  PHActsGSF(const std::string& name = "PHActsGSF");

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_pp_mode(bool mode) {m_pp_mode = mode;}

  void useActsEvaluator(bool actsEvaluator)
  {
    m_actsEvaluator = actsEvaluator;
  }

 private:
  int getNodes(PHCompositeNode* topNode);
  std::shared_ptr<Acts::PerigeeSurface> makePerigee(SvtxTrack* track) const;
  ActsTrackFittingAlgorithm::TrackParameters makeSeed(
      SvtxTrack* track,
      const std::shared_ptr<Acts::PerigeeSurface>& psurf) const;
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
                   ActsTrackFittingAlgorithm::TrackContainer& tracks,
                   const TrackSeed* seed, const ActsTrackFittingAlgorithm::MeasurementContainer& measurements);
  void updateSvtxTrack(std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
                       Trajectory::IndexedParameters& paramsMap,
                       ActsTrackFittingAlgorithm::TrackContainer& tracks,
                       SvtxTrack* track);
  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track);

  ActsGeometry* m_tGeometry = nullptr;
  TrkrClusterContainer* m_clusterContainer = nullptr;
  SvtxTrackMap* m_trackMap = nullptr;
  SvtxVertexMap* m_vertexMap = nullptr;

//  alignmentTransformationContainer* m_alignmentTransformationMap = nullptr;  // added for testing purposes
  alignmentTransformationContainer* m_alignmentTransformationMapTransient = nullptr;
  std::set< Acts::GeometryIdentifier> m_transient_id_set;
  Acts::GeometryContext m_transient_geocontext;

  // Tpc Global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  std::string m_trackMapName = "SvtxTrackMap";
  std::string _seed_track_map_name = "SeedTrackMap";
//  unsigned int m_pHypothesis = 11;

  bool m_pp_mode = false;

  bool m_actsEvaluator = false;
  std::unique_ptr<ActsEvaluator> m_evaluator = nullptr;
  std::string m_evalname = "ActsEvaluator.root";

  SvtxTrackMap* m_seedTracks = nullptr;

  ClusterErrorPara _ClusErrPara;

  ActsTrackFittingAlgorithm::Config m_fitCfg;
};

#endif  // PHACTSGSF_H
