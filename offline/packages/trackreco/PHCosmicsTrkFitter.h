
#ifndef TRACKRECO_PHCOSMICSTRKFITTER_H
#define TRACKRECO_PHCOSMICSTRKFITTER_H

#include "ActsAlignmentStates.h"
#include "ActsEvaluator.h"

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSourceLink.h>
#include <trackbase/ActsTrackFittingAlgorithm.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/alignmentTransformationContainer.h>

#include <tpc/TpcClusterMover.h>
#include <tpc/TpcClusterZCrossingCorrection.h>
#include <tpc/TpcDistortionCorrection.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <ActsExamples/EventData/Trajectories.hpp>

#include <memory>
#include <string>

#include <TFile.h>
#include <TTree.h>

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

class PHCosmicsTrkFitter : public SubsysReco
{
 public:
  /// Default constructor
  PHCosmicsTrkFitter(const std::string& name = "PHCosmicsTrkFitter");

  /// Destructor
  ~PHCosmicsTrkFitter() override = default;

  /// End, write and close files
  int End(PHCompositeNode* topNode) override;

  /// Get and create nodes
  int InitRun(PHCompositeNode* topNode) override;

  /// Process each event by calling the fitter
  int process_event(PHCompositeNode* topNode) override;

  int ResetEvent(PHCompositeNode* topNode) override;

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
  void seedAnalysis() { m_seedClusAnalysis = true; }
  void commissioning(bool com) { m_commissioning = com; }

  void useOutlierFinder(bool outlier) { m_useOutlierFinder = outlier; }

  void SetIteration(int iter) { _n_iteration = iter; }
  void set_track_map_name(const std::string& map_name) { _track_map_name = map_name; }
  void set_seed_track_map_name(const std::string& map_name) { _seed_track_map_name = map_name; }

  void ignoreLayer(int layer) { m_ignoreLayer.insert(layer); }
  void setVertexRadius(const float rad) { m_vertexRadius = rad; }

 private:
  /// Get all the nodes
  int getNodes(PHCompositeNode* topNode);

  /// Create new nodes
  int createNodes(PHCompositeNode* topNode);

  void loopTracks(Acts::Logging::Level logLevel);
  SourceLinkVec getSourceLinks(TrackSeed* track,
                               ActsTrackFittingAlgorithm::MeasurementContainer& measurements,
                               short int crossing,
                               int& charge,
                               float& cosmicslope);

  /// Convert the acts track fit result to an svtx track
  void updateSvtxTrack(std::vector<Acts::MultiTrajectoryTraits::IndexType>& tips,
                       Trajectory::IndexedParameters& paramsMap,
                       ActsTrackFittingAlgorithm::TrackContainer& tracks,
                       SvtxTrack* track);

  /// Helper function to call either the regular navigation or direct
  /// navigation, depending on m_fitSiliconMMs
  inline ActsTrackFittingAlgorithm::TrackFitterResult fitTrack(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const ActsTrackFittingAlgorithm::TrackParameters& seed,
      const ActsTrackFittingAlgorithm::GeneralFitterOptions&
          kfOptions,
      const CalibratorAdapter& calibrator,
      ActsTrackFittingAlgorithm::TrackContainer& tracks);

  bool getTrackFitResult(FitResult& fitOutput, TrackSeed* seed,
                         SvtxTrack* track,
                         ActsTrackFittingAlgorithm::TrackContainer& tracks,
                         const ActsTrackFittingAlgorithm::MeasurementContainer& measurements);

  Acts::BoundSquareMatrix setDefaultCovariance() const;
  void printTrackSeed(const ActsTrackFittingAlgorithm::TrackParameters& seed) const;
  void makeBranches();

  /// Event counter
  int m_event = 0;

  /// Options that Acts::Fitter needs to run from MakeActsGeometry
  ActsGeometry* m_tGeometry = nullptr;

  /// Configuration containing the fitting function instance
  ActsTrackFittingAlgorithm::Config m_fitCfg;

  /// TrackMap containing SvtxTracks
  alignmentTransformationContainer* m_alignmentTransformationMap = nullptr;  // added for testing purposes
  SvtxTrackMap* m_trackMap = nullptr;
  SvtxTrackMap* m_directedTrackMap = nullptr;
  TrkrClusterContainer* m_clusterContainer = nullptr;
  TrackSeedContainer* m_seedMap = nullptr;
  TrackSeedContainer* m_tpcSeeds = nullptr;
  TrackSeedContainer* m_siliconSeeds = nullptr;

  /// Number of acts fits that returned an error
  int m_nBadFits = 0;

  //! bool to fill alignment state map for further processing
  bool m_commissioning = true;

  /// A bool to update the SvtxTrackState information (or not)
  bool m_fillSvtxTrackStates = true;

  /// A bool to use the chi2 outlier finder in the track fitting
  bool m_useOutlierFinder = false;
  ResidualOutlierFinder m_outlierFinder;

  float m_vertexRadius = 80;

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

  // cluster mover utility class
  TpcClusterMover _clusterMover;
  ClusterErrorPara _ClusErrPara;

  std::set<int> m_ignoreLayer;

  std::string m_fieldMap = "";


  int _n_iteration = 0;
  std::string _track_map_name = "SvtxTrackMap";
  std::string _seed_track_map_name = "SeedTrackMap";

  /// Default particle assumption to muon
  unsigned int m_pHypothesis = 13;

  SvtxAlignmentStateMap* m_alignmentStateMap = nullptr;
  ActsAlignmentStates m_alignStates;


  //! for diagnosing seed param + clusters
  bool m_seedClusAnalysis = false;
  std::unique_ptr<TFile> m_outfile = nullptr;
  std::unique_ptr<TTree> m_tree = nullptr;
  int m_seed = std::numeric_limits<int>::max();
  float m_R = NAN;
  float m_X0 = NAN;
  float m_Y0 = NAN;
  float m_Z0 = NAN;
  float m_slope = NAN;
  float m_pcax = NAN;
  float m_pcay = NAN;
  float m_pcaz = NAN;
  float m_px = NAN;
  float m_py = NAN;
  float m_pz = NAN;
  int m_charge = std::numeric_limits<int>::max();
  int m_nmaps = std::numeric_limits<int>::max();
  int m_nintt = std::numeric_limits<int>::max();
  int m_ntpc = std::numeric_limits<int>::max();
  int m_nmm = std::numeric_limits<int>::max();
  std::vector<float> m_locx, m_locy, m_x, m_y, m_z, m_r, m_layer,m_phi, m_eta, 
    m_phisize, m_zsize, m_ephi, m_ez;
  void clearVectors();
  void fillVectors(TrackSeed* tpcseed, TrackSeed *siseed);
  ClusterErrorPara m_clusErrPara;
};

#endif
