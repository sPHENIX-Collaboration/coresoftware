/*!
 *  \file		PHActsTrkFitter.h
 *  \brief		Refit SvtxTracks with Acts.
 *  \details	Refit SvtxTracks with Acts
 *  \author		Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_ACTSTRKFITTER_H
#define TRACKRECO_ACTSTRKFITTER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsSurfaceMaps.h>

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ActsExamples/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

#include <boost/bimap.hpp>

#include <memory>
#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

namespace ActsExamples
{
  class TrkrClusterSourceLink;
}

class MakeActsGeometry;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class TrkrClusterContainer;

using SourceLink = ActsExamples::TrkrClusterSourceLink;
using FitResult = Acts::KalmanFitterResult<SourceLink>;
using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;
using Measurement = Acts::Measurement<ActsExamples::TrkrClusterSourceLink,
                                      Acts::BoundIndices,
                                      Acts::eBoundLoc0,
                                      Acts::eBoundLoc1>;
using SurfacePtrVec = std::vector<const Acts::Surface*>;
using SourceLinkVec = std::vector<SourceLink>;

class PHActsTrkFitter : public SubsysReco
{
 public:
  /// Default constructor
  PHActsTrkFitter(const std::string& name = "PHActsTrkFitter");

  /// Destructor
  ~PHActsTrkFitter() override;

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

  void setUpdateSvtxTrackStates(bool fillSvtxTrackStates)
       { m_fillSvtxTrackStates = fillSvtxTrackStates; }   

 private:

  /// Event counter
  int m_event;

  /// Get all the nodes
  int getNodes(PHCompositeNode *topNode);

  /// Create new nodes
  int createNodes(PHCompositeNode *topNode);

  void loopTracks(Acts::Logging::Level logLevel);
  SourceLinkVec getSourceLinks(SvtxTrack *track);
  Acts::Vector3D getVertex(SvtxTrack *track);

  /// Convert the acts track fit result to an svtx track
  void updateSvtxTrack(Trajectory traj, SvtxTrack* track);

  /// Helper function to call either the regular navigation or direct
  /// navigation, depending on m_fitSiliconMMs
  ActsExamples::TrkrClusterFittingAlgorithm::FitterResult fitTrack(
           const SourceLinkVec& sourceLinks, 
	   const ActsExamples::TrackParameters& seed,
	   const Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>& 
	         kfOptions,
	   const SurfacePtrVec& surfSequence);

  /// Functions to get list of sorted surfaces for direct navigation, if
  /// applicable
  SourceLinkVec getSurfaceVector(SourceLinkVec sourceLinks, 
				 SurfacePtrVec& surfaces);
  void checkSurfaceVec(SurfacePtrVec& surfaces);
  void getTrackFitResult(const FitResult& fitOutput, 
			 SvtxTrack* track);

  Surface getSurface(TrkrDefs::cluskey cluskey,TrkrDefs::subsurfkey surfkey);
  Surface getSiliconSurface(TrkrDefs::hitsetkey hitsetkey);
  Surface getTpcSurface(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::subsurfkey surfkey);
  Surface getMMSurface(TrkrDefs::hitsetkey hitsetkey);

  Acts::BoundSymMatrix setDefaultCovariance();
  void printTrackSeed(ActsExamples::TrackParameters seed);

  /// Options that Acts::Fitter needs to run from MakeActsGeometry
  ActsTrackingGeometry *m_tGeometry;

  /// Configuration containing the fitting function instance
  ActsExamples::TrkrClusterFittingAlgorithm::Config m_fitCfg;

  /// TrackMap containing SvtxTracks
  SvtxTrackMap *m_trackMap, *m_directedTrackMap;
  SvtxVertexMap *m_vertexMap;
  TrkrClusterContainer *m_clusterContainer;
  ActsSurfaceMaps *m_surfMaps;

  /// Number of acts fits that returned an error
  int m_nBadFits;

  /// Boolean to use normal tracking geometry navigator or the
  /// Acts::DirectedNavigator with a list of sorted silicon+MM surfaces
  bool m_fitSiliconMMs;

  /// A bool to update the SvtxTrackState information (or not)
  bool m_fillSvtxTrackStates;

  /// Variables for doing event time execution analysis
  bool m_timeAnalysis;
  TFile *m_timeFile;
  TH1 *h_eventTime;
  TH2 *h_fitTime;
  TH1 *h_updateTime;
  TH1 *h_stateTime;
  TH1 *h_rotTime;
};

#endif
