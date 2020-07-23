/*!
 *  \file		PHActsTrkFitter.h
 *  \brief		Refit SvtxTracks with Acts.
 *  \details	Refit SvtxTracks with Acts
 *  \author		Tony Frawley <afrawley@fsu.edu>
 */

#ifndef TRACKRECO_ACTSTRKFITTER_H
#define TRACKRECO_ACTSTRKFITTER_H

#include "PHTrackFitting.h"
#include "PHActsSourceLinks.h"

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ACTFW/Fitting/TrkrClusterFittingAlgorithm.hpp>
#include <ACTFW/EventData/TrkrClusterMultiTrajectory.hpp>

#include <memory>
#include <string>
#include <TFile.h>
#include <TH1.h>

namespace FW
{
  namespace Data
  {
    class TrkrClusterSourceLink;
  }
}  // namespace FW

class ActsTrack;
class MakeActsGeometry;
class SvtxTrack;
class SvtxTrackMap;

using SourceLink = FW::Data::TrkrClusterSourceLink;
using FitResult = Acts::KalmanFitterResult<SourceLink>;
using Trajectory = FW::TrkrClusterMultiTrajectory;
using Measurement = Acts::Measurement<FW::Data::TrkrClusterSourceLink,
                                      Acts::ParDef::eLOC_0,
                                      Acts::ParDef::eLOC_1>;
class PHActsTrkFitter : public PHTrackFitting
{
 public:
  /// Default constructor
  PHActsTrkFitter(const std::string& name = "PHActsTrkFitter");

  /// Destructor
  ~PHActsTrkFitter();

  /// End, write and close files
  int End(PHCompositeNode*);

  /// Get and create nodes
  int Setup(PHCompositeNode* topNode);

  /// Process each event by calling the fitter
  int Process();

  int ResetEvent(PHCompositeNode *topNode);

  void setTimeAnalysis(bool time){m_timeAnalysis = time;}

 private:

  // DCA 
  void calculateDCA(const Acts::BoundParameters param);

  void fillSvtxTrackStates(const Trajectory traj, const size_t &trackTip, SvtxTrack *svtx_track);

  TrkrDefs::cluskey getClusKey(const unsigned int hitID);

  /// Event counter
  int m_event;

  /// Get all the nodes
  int getNodes(PHCompositeNode*);

  /// Create new nodes
  int createNodes(PHCompositeNode*);

  /// Convert the acts track fit result to an svtx track
  void updateSvtxTrack(Trajectory traj, const unsigned int trackKey);

  /// Map of Acts fit results and track key to be placed on node tree
  std::map<const unsigned int, Trajectory> 
    *m_actsFitResults;

  /// Map of acts tracks and track key created by PHActsTracks
  std::map<unsigned int, ActsTrack>* m_actsProtoTracks;

  /// Options that Acts::Fitter needs to run from MakeActsGeometry
  ActsTrackingGeometry *m_tGeometry;

  /// Configuration containing the fitting function instance
  FW::TrkrClusterFittingAlgorithm::Config fitCfg;

  /// TrackMap containing SvtxTracks
  SvtxTrackMap *m_trackMap;

  // map relating acts hitid's to clusterkeys
  std::map<TrkrDefs::cluskey, unsigned int> *m_hitIdClusKey{nullptr};

  int m_nBadFits;

  float m_dca3Dxy;          /// fitted parameter 3D DCA in xy plane
  float m_dca3Dz;           /// fitted parameter 3D DCA in z plane
  float m_dca3DxyCov;       /// fitted parameter 3D DCA covariance in xy
  float m_dca3DzCov;        /// fitted parameter 3D DCA covariance in z


  /// Variables for doing event time execution analysis
  bool m_timeAnalysis;
  TFile *m_timeFile;
  TH1 *h_eventTime;
  
};

#endif
