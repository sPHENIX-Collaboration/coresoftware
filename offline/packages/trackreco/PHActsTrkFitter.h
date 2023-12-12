/*!
 *  \file		PHActsTrkFitter.h
 *  \brief		Refit SvtxTracks with Acts.
 *  \details	        Refit SvtxTracks with Acts
 *  \author		Joe Osborn, Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_ACTSTRKFITTER_H
#define TRACKRECO_ACTSTRKFITTER_H

#include "ActsAlignmentStates.h"
#include "ActsEvaluator.h"

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSourceLink.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/alignmentTransformationContainer.h>

#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <memory>
#include <string>

#include <trackbase/alignmentTransformationContainer.h>

class MakeActsGeometry;
class SvtxTrack;
class SvtxTrackMap;
class TrackSeed;
class TrackSeedContainer;
class TrkrClusterContainer;
class TpcDistortionCorrectionContainer;
class SvtxAlignmentStateMap;

using SourceLink = ActsSourceLink;
using FitResult = ActsTrackFittingAlgorithm::TrackFitterResult;
using Trajectory = ActsExamples::Trajectories;
using Measurement = Acts::Measurement<Acts::BoundIndices, 2>;
using SurfacePtrVec = std::vector<const Acts::Surface*>;
using SourceLinkVec = std::vector<Acts::SourceLink>;

class PHActsTrkFitter : public SubsysReco
{
 public:
  /// Default constructor
  PHActsTrkFitter(const std::string& name = "PHActsTrkFitter");

  /// Destructor
  ~PHActsTrkFitter() override = default;

  /// End, write and close files
  int End(PHCompositeNode* topNode) override;

  /// Get and create nodes
  int InitRun(PHCompositeNode* topNode) override;

  /// Process each event by calling the fitter
  int process_event(PHCompositeNode* topNode) override;

  int ResetEvent(PHCompositeNode* topNode) override;

  /// Do some internal time benchmarking analysis
  void doTimeAnalysis(bool timeAnalysis) { m_timeAnalysis = timeAnalysis; }

  /// Run the direct navigator to fit only tracks with silicon+MM hits
  void fitSiliconMMs(bool fitSiliconMMs)
  {
    m_fitSiliconMMs = fitSiliconMMs;
  }

  /// require micromegas in SiliconMM fits
  void setUseMicromegas(bool value)
  {
    m_useMicromegas = value;
  }

  void setUpdateSvtxTrackStates(bool fillSvtxTrackStates)
  {
    m_fillSvtxTrackStates = fillSvtxTrackStates;
  }

  void useActsEvaluator(bool actsEvaluator)
  {
    m_actsEvaluator = actsEvaluator;
  }

  void setEvaluatorName(const std::string &name) { m_evalname = name; }
  void setFieldMap(const std::string &fieldMap)
  {
    m_fieldMap = fieldMap;
  }

  void setAbsPdgHypothesis(unsigned int pHypothesis)
  {
    m_pHypothesis = pHypothesis;
  }

  void commissioning(bool com) { m_commissioning = com; }

  void useOutlierFinder(bool outlier) { m_useOutlierFinder = outlier; }

  void SetIteration(int iter) { _n_iteration = iter; }
  void set_track_map_name(const std::string& map_name) { _track_map_name = map_name; }
  void set_seed_track_map_name(const std::string& map_name) { _seed_track_map_name = map_name; }

  /// Set flag for pp running
  void set_pp_mode(bool ispp) { m_pp_mode = ispp; }

  void ignoreLayer(int layer) { m_ignoreLayer.insert(layer); }

 private:
  /// Get all the nodes
  int getNodes(PHCompositeNode* topNode);

  /// Create new nodes
  int createNodes(PHCompositeNode* topNode);

  void loopTracks(Acts::Logging::Level logLevel);
  SourceLinkVec getSourceLinks(TrackSeed* track,
                               ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
                               short int crossing);

  /// Convert the acts track fit result to an svtx track
  void updateSvtxTrack(std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
                       Trajectory::IndexedParameters& paramsMap,
                       ActsTrackFittingAlgorithm::TrackContainer& tracks,
                       SvtxTrack* track);

  /// Helper function to call either the regular navigation or direct
  /// navigation, depending on m_fitSiliconMMs
  ActsTrackFittingAlgorithm::TrackFitterResult fitTrack(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const ActsTrackFittingAlgorithm::TrackParameters& seed,
      const ActsTrackFittingAlgorithm::GeneralFitterOptions&
          kfOptions,
      const SurfacePtrVec& surfSequence,
      const CalibratorAdapter& calibrator,
      ActsTrackFittingAlgorithm::TrackContainer& tracks);

  /// Functions to get list of sorted surfaces for direct navigation, if
  /// applicable
  SourceLinkVec getSurfaceVector(const SourceLinkVec& sourceLinks,
                                 SurfacePtrVec& surfaces) const;
  void checkSurfaceVec(SurfacePtrVec& surfaces) const;

  bool getTrackFitResult(FitResult& fitOutput, TrackSeed* seed,
                         SvtxTrack* track,
                         ActsTrackFittingAlgorithm::TrackContainer& tracks,
                         const ActsTrackFittingAlgorithm::MeasurementContainer& measurements);

  Acts::BoundSquareMatrix setDefaultCovariance() const;
  void printTrackSeed(const ActsTrackFittingAlgorithm::TrackParameters& seed) const;

  /// Event counter
  int m_event = 0;

  /// Options that Acts::Fitter needs to run from MakeActsGeometry
  ActsGeometry* m_tGeometry = nullptr;

  /// Configuration containing the fitting function instance
  ActsTrackFittingAlgorithm::Config m_fitCfg;

  /// TrackMap containing SvtxTracks
  alignmentTransformationContainer* m_alignmentTransformationMap = nullptr;  // added for testing purposes
  alignmentTransformationContainer* m_alignmentTransformationMapTransient = nullptr;  
  std::set< Acts::GeometryIdentifier> m_transient_id_set;
  Acts::GeometryContext m_transient_geocontext;
  SvtxTrackMap* m_trackMap = nullptr;
  SvtxTrackMap* m_directedTrackMap = nullptr;
  TrkrClusterContainer* m_clusterContainer = nullptr;
  TrackSeedContainer* m_seedMap = nullptr;
  TrackSeedContainer* m_tpcSeeds = nullptr;
  TrackSeedContainer* m_siliconSeeds = nullptr;

  /// Number of acts fits that returned an error
  int m_nBadFits = 0;

  /// Boolean to use normal tracking geometry navigator or the
  /// Acts::DirectedNavigator with a list of sorted silicon+MM surfaces
  bool m_fitSiliconMMs = false;

  /// requires micromegas present when fitting silicon-MM surfaces
  bool m_useMicromegas = true;

  /// A bool to update the SvtxTrackState information (or not)
  bool m_fillSvtxTrackStates = true;

  /// A bool to use the chi2 outlier finder in the track fitting
  bool m_useOutlierFinder = false;
  ResidualOutlierFinder m_outlierFinder;

  /// Flag for pp running
  bool m_pp_mode = false;

  // max variation of bunch crossing away from crossing_estimate
  short int max_bunch_search = 2;

  bool m_actsEvaluator = false;
  std::unique_ptr<ActsEvaluator> m_evaluator = nullptr;
  std::string m_evalname = "ActsEvaluator.root";

  std::map<const unsigned int, Trajectory>* m_trajectories = nullptr;
  SvtxTrackMap* m_seedTracks = nullptr;

  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};

  /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  ClusterErrorPara _ClusErrPara;

  std::set<int> m_ignoreLayer;

  std::string m_fieldMap = "";

  int _n_iteration = 0;
  std::string _track_map_name = "SvtxTrackMap";
  std::string _seed_track_map_name = "SeedTrackMap";

  /// Default particle assumption to pion
  unsigned int m_pHypothesis = 211;

  SvtxAlignmentStateMap* m_alignmentStateMap = nullptr;
  ActsAlignmentStates m_alignStates;
  bool m_commissioning = false;

  /// Variables for doing event time execution analysis
  bool m_timeAnalysis = false;
  TFile* m_timeFile = nullptr;
  TH1* h_eventTime = nullptr;
  TH2* h_fitTime = nullptr;
  TH1* h_updateTime = nullptr;
  TH1* h_stateTime = nullptr;
  TH1* h_rotTime = nullptr;
};

#endif
