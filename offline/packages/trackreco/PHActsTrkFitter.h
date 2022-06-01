/*!
 *  \file		PHActsTrkFitter.h
 *  \brief		Refit SvtxTracks with Acts.
 *  \details	Refit SvtxTracks with Acts
 *  \author		Joe Osborn, Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_ACTSTRKFITTER_H
#define TRACKRECO_ACTSTRKFITTER_H


#include <fun4all/SubsysReco.h>

#include "ResidualOutlierFinder.h"

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

#include <memory>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


class MakeActsGeometry;
class SvtxTrack;
class SvtxTrack_v4;
class SvtxTrackMap;
class TrackSeed;
class TrackSeedContainer;
class TrkrClusterContainer;
class TrkrClusterIterationMapv1;

using SourceLink = ActsExamples::IndexSourceLink;
using FitResult = Acts::KalmanFitterResult;
using Trajectory = ActsExamples::Trajectories;
using Measurement = Acts::Measurement<Acts::BoundIndices,2>;
using SurfacePtrVec = std::vector<const Acts::Surface*>;
using SourceLinkVec = std::vector<SourceLink>;

class PHActsTrkFitter : public SubsysReco
{
 public:
  /// Default constructor
  PHActsTrkFitter(const std::string& name = "PHActsTrkFitter");

  /// Destructor
  ~PHActsTrkFitter() override = default;

  /// End, write and close files
  int End(PHCompositeNode *topNode) override;

  /// Get and create nodes
  int InitRun(PHCompositeNode* topNode) override;

  /// Process each event by calling the fitter
  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  /// Do some internal time benchmarking analysis
  void doTimeAnalysis(bool timeAnalysis){m_timeAnalysis = timeAnalysis;}

  /// Run the direct navigator to fit only tracks with silicon+MM hits
  void fitSiliconMMs(bool fitSiliconMMs)
       {m_fitSiliconMMs = fitSiliconMMs;}

  /// require micromegas in SiliconMM fits
  void setUseMicromegas( bool value )
  { m_useMicromegas = value; }

  void setUpdateSvtxTrackStates(bool fillSvtxTrackStates)
       { m_fillSvtxTrackStates = fillSvtxTrackStates; }   

  void useActsEvaluator(bool actsEvaluator)
  { m_actsEvaluator = actsEvaluator; }
  
  void setFieldMap(std::string& fieldMap)
  { m_fieldMap = fieldMap; }

  void setAbsPdgHypothesis(unsigned int pHypothesis)
  { m_pHypothesis = pHypothesis; }

  void useOutlierFinder(bool outlier) { m_useOutlierFinder = outlier; }

  void SetIteration(int iter){_n_iteration = iter;}
  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void set_seed_track_map_name(const std::string &map_name) { _seed_track_map_name = map_name; }

 private:

  /// Get all the nodes
  int getNodes(PHCompositeNode *topNode);

  /// Create new nodes
  int createNodes(PHCompositeNode *topNode);

  void loopTracks(Acts::Logging::Level logLevel);
  SourceLinkVec getSourceLinks(TrackSeed *track, 
			       ActsExamples::MeasurementContainer& measurements,
			       short int crossing);

  /// Convert the acts track fit result to an svtx track
  void updateSvtxTrack(Trajectory traj, std::unique_ptr<SvtxTrack_v4>& track);

  /// Helper function to call either the regular navigation or direct
  /// navigation, depending on m_fitSiliconMMs
  ActsExamples::TrackFittingAlgorithm::TrackFitterResult fitTrack(
           const std::vector<std::reference_wrapper<const SourceLink>>& sourceLinks, 
	   const ActsExamples::TrackParameters& seed,
	   const ActsExamples::TrackFittingAlgorithm::TrackFitterOptions& 
	     kfOptions,
	   const SurfacePtrVec& surfSequence);

  /// Functions to get list of sorted surfaces for direct navigation, if
  /// applicable
  SourceLinkVec getSurfaceVector(const SourceLinkVec& sourceLinks, 
				 SurfacePtrVec& surfaces) const;
  void checkSurfaceVec(SurfacePtrVec& surfaces) const;

  bool getTrackFitResult(const FitResult& fitOutput, std::unique_ptr<SvtxTrack_v4>& track);

  Surface getSurface(TrkrDefs::cluskey cluskey,TrkrDefs::subsurfkey surfkey) const;
  Surface getSiliconSurface(TrkrDefs::hitsetkey hitsetkey) const;
  Surface getTpcSurface(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::subsurfkey surfkey) const;
  Surface getMMSurface(TrkrDefs::hitsetkey hitsetkey) const;

  Surface get_tpc_surface_from_coords(TrkrDefs::hitsetkey hitsetkey,
						       Acts::Vector3 world,
						       ActsSurfaceMaps *surfMaps,
						       ActsTrackingGeometry *tGeometry,
						       TrkrDefs::subsurfkey& subsurfkey);

    Acts::BoundSymMatrix setDefaultCovariance() const;
  void printTrackSeed(const ActsExamples::TrackParameters& seed) const;

  /// Event counter
  int m_event = 0;

  /// Options that Acts::Fitter needs to run from MakeActsGeometry
  ActsTrackingGeometry *m_tGeometry = nullptr;

  /// Configuration containing the fitting function instance
  ActsExamples::TrackFittingAlgorithm::Config m_fitCfg;

  /// TrackMap containing SvtxTracks
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxTrackMap *m_directedTrackMap = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  ActsSurfaceMaps *m_surfMaps = nullptr;
  TrackSeedContainer *m_seedMap = nullptr;
  TrackSeedContainer *m_tpcSeeds = nullptr;
  TrackSeedContainer *m_siliconSeeds = nullptr;
  
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

  bool m_actsEvaluator = false;
  std::map<const unsigned int, Trajectory> *m_trajectories = nullptr;
  SvtxTrackMap *m_seedTracks = nullptr;

  TpcClusterZCrossingCorrection m_clusterCrossingCorrection;
  TpcDistortionCorrectionContainer* _dcc_static{nullptr};
  TpcDistortionCorrectionContainer* _dcc_average{nullptr};
  TpcDistortionCorrectionContainer* _dcc_fluctuation{nullptr};

 /// tpc distortion correction utility class
  TpcDistortionCorrection _distortionCorrection;

  // cluster mover utility class
  TpcClusterMover _clusterMover;

  std::string m_fieldMap = "";
  TrkrClusterIterationMapv1* _iteration_map = nullptr;
  int _n_iteration = 0;
  std::string _track_map_name = "SvtxTrackMap";
  std::string _seed_track_map_name = "SeedTrackMap";

  /// Default particle assumption to pion
  unsigned int m_pHypothesis = 211;

  /// Variables for doing event time execution analysis
  bool m_timeAnalysis = false;
  TFile *m_timeFile = nullptr;
  TH1 *h_eventTime = nullptr;
  TH2 *h_fitTime = nullptr;
  TH1 *h_updateTime = nullptr;
  TH1 *h_stateTime = nullptr;
  TH1 *h_rotTime = nullptr;
};

#endif
