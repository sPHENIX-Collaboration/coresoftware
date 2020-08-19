#ifndef TRACKRECO_PHACTSTRKPROP_H
#define TRACKRECO_PHACTSTRKPROP_H

#include "PHTrackPropagating.h"
#include "PHActsSourceLinks.h"
#include "ActsTrack.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/Geometry/GeometryID.hpp>

#include <Acts/TrackFinder/CKFSourceLinkSelector.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>

#include <ACTFW/EventData/TrkrClusterSourceLink.hpp>
#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/EventData/TrkrClusterMultiTrajectory.hpp>

#include <ACTFW/TrackFinding/TrkrClusterFindingAlgorithm.hpp>

#include <memory>
#include <string>
#include <map>

class PHCompositeNode;
class ActsTrack;
class MakeActsGeometry;
class SvtxTrack;
class SvtxTrackMap;

class TFile;
class TH1;

namespace FW
{
  namespace Data
  {
    class TrkrClusterSourceLink;
  }
}
namespace Acts
{
  class Surface;
  class PerigeeSurface;
}

using PerigeeSurface = std::shared_ptr<const Acts::PerigeeSurface>;
using Surface = std::shared_ptr<const Acts::Surface>;
using SourceLink = FW::Data::TrkrClusterSourceLink;

using SourceLinkSelector = Acts::CKFSourceLinkSelector;
using SourceLinkSelectorConfig = typename SourceLinkSelector::Config;

using CKFFitResult = Acts::CombinatorialKalmanFilterResult<SourceLink>;
using Trajectory = FW::TrkrClusterMultiTrajectory;

using Measurement = Acts::Measurement<FW::Data::TrkrClusterSourceLink,
                                      Acts::BoundParametersIndices,
                                      Acts::ParDef::eLOC_0,
                                      Acts::ParDef::eLOC_1>;

class PHActsTrkProp : public PHTrackPropagating
{
 public:
  /// Default constructor
  PHActsTrkProp(const std::string& name = "PHActsTrkProp");

  /// Destructor
  ~PHActsTrkProp();

  /// End, write and close files
  int End();

  /// Get and create nodes
  int Setup(PHCompositeNode* topNode);

  /// Process each event by calling the fitter
  int Process();

  /// Reset maps event by event
  int ResetEvent(PHCompositeNode *topNode);

  void doTimeAnalysis(bool timeAnalysis) { m_timeAnalysis = timeAnalysis;}

 private:
  /// Event counter
  int m_event;

  bool m_timeAnalysis;
  
  TFile *m_timeFile;
  TH1 *h_eventTime;

  /// Num bad fit counter
  int m_nBadFits;

  /// Get all the nodes
  int getNodes(PHCompositeNode *topNode);

  /// Create new nodes
  void createNodes(PHCompositeNode *topNode);

  /// Helper function to make an Acts::GeometryID for SL selection
  Acts::GeometryID makeId(int volume = 0, 
			  int layer = 0, 
			  int sensitive = 0);

  /// Wipe and recreate the SvtxTrackMap with Acts output
  void updateSvtxTrack(Trajectory traj, 
		       const unsigned int trackKey, 
		       Acts::Vector3D vertex);

  /// Get all source links in a given event
  std::vector<SourceLink> getEventSourceLinks();

  ActsTrackingGeometry *m_tGeometry;

  /// Track map with Svtx objects
  SvtxTrackMap *m_trackMap;

  /// Track map with ActsTracks, created in PHActsTracks
  std::map<unsigned int, ActsTrack> *m_actsProtoTracks;
  
  /// Acts MultiTrajectories for ActsEvaluator
  std::map<const unsigned int, Trajectory> *m_actsFitResults;

  /// Map that correlates track key with track tip for ActsEvaluator
  /// Identifiers are <TrajNum, <trackTip, SvtxTrackKey>> where 
  /// TrajNum has a one-to-one map to the original seed SvtxTrackKey
  std::map<const unsigned int, 
    std::map<const size_t, const unsigned int>> *m_actsTrackKeyMap;

  /// Map of cluster keys to hit ids, for identifying clusters belonging to track
  std::map<TrkrDefs::cluskey, unsigned int> *m_hitIdClusKey;

  /// Acts source links created by PHActsSourceLinks
  /// SourceLink is defined as TrkrClusterSourceLink elsewhere
  std::map<unsigned int, SourceLink> *m_sourceLinks;

  PHCompositeNode *m_topNode;

  /// SourceLinkSelector to help the CKF identify which source links 
  /// may belong to a track seed based on geometry considerations
  SourceLinkSelectorConfig m_sourceLinkSelectorConfig;
 
  /// Configuration containing the finding function instance
  FW::TrkrClusterFindingAlgorithm::Config findCfg;


};

#endif
