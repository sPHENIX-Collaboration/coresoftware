#ifndef TRACKRECO_PHACTSTRKPROP_H
#define TRACKRECO_PHACTSTRKPROP_H

#include "PHTrackPropagating.h"
#include "ActsTrackingGeometry.h"
#include "ActsTrack.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/Geometry/GeometryIdentifier.hpp>

#include <Acts/TrackFinding/CKFSourceLinkSelector.hpp>

#include <Acts/EventData/MeasurementHelpers.hpp>

#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

#include <ActsExamples/TrackFinding/TrkrClusterFindingAlgorithm.hpp>

#include <boost/bimap.hpp>

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

namespace ActsExamples
{
  class TrkrClusterSourceLink;  
}
namespace Acts
{
  class Surface;
  class PerigeeSurface;
}

using PerigeeSurface = std::shared_ptr<const Acts::PerigeeSurface>;
using Surface = std::shared_ptr<const Acts::Surface>;
using SourceLink = ActsExamples::TrkrClusterSourceLink;

using SourceLinkSelector = Acts::CKFSourceLinkSelector;
using SourceLinkSelectorConfig = typename SourceLinkSelector::Config;

using CKFFitResult = Acts::CombinatorialKalmanFilterResult<SourceLink>;
using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;

using Measurement = Acts::Measurement<ActsExamples::TrkrClusterSourceLink,
                                      Acts::BoundIndices,
                                      Acts::eBoundLoc0,
                                      Acts::eBoundLoc1>;

typedef boost::bimap<TrkrDefs::cluskey, unsigned int> CluskeyBimap;

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

  void setVolumeMaxChi2(const int vol, const float maxChi2);
  void setVolumeLayerMaxChi2(const int vol, const int layer,
			     const float maxChi2);
  void resetCovariance(bool resetCovariance){m_resetCovariance = resetCovariance;}

 private:

  /// Event counter
  int m_event;

  bool m_timeAnalysis; 
  TFile *m_timeFile;
  TH1 *h_eventTime;

  void initializeLayerSelector();

  /// Map to hold maximum allowable measurement chi2 in each
  /// volume identifier in Acts
  std::map<const int, const float> m_volMaxChi2;

  /// array of maps to hold maximum allowable measurement chi 2
  /// in a particular layer of a particular volume id in acts
  /// Entries are 0 - MVTX, 1 - INTT, 2 - TPC
  std::vector<std::map<const int, const float>> m_volLayerMaxChi2;

  /// Num bad fit counter
  int m_nBadFits;

  /// Get all the nodes
  int getNodes(PHCompositeNode *topNode);

  /// Create new nodes
  void createNodes(PHCompositeNode *topNode);

  /// Helper function to make an Acts::GeometryID for SL selection
  Acts::GeometryIdentifier makeId(int volume = 0, 
				  int layer = 0, 
				  int sensitive = 0);

  /// Wipe and recreate the SvtxTrackMap with Acts output
  void updateSvtxTrack(Trajectory traj, 
		       const unsigned int trackKey, 
		       Acts::Vector3D vertex);

  /// Get all source links in a given event
  std::vector<SourceLink> getEventSourceLinks();

  /// Setup the source link selector criteria
  void setupSourceLinkSelection();

  ActsTrackingGeometry *m_tGeometry;

  bool m_resetCovariance;

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
  CluskeyBimap *m_hitIdClusKey;

  /// Acts source links created by PHActsSourceLinks
  /// SourceLink is defined as TrkrClusterSourceLink elsewhere
  std::map<unsigned int, SourceLink> *m_sourceLinks;

  PHCompositeNode *m_topNode;

  /// SourceLinkSelector to help the CKF identify which source links 
  /// may belong to a track seed based on geometry considerations
  SourceLinkSelectorConfig m_sourceLinkSelectorConfig;
 
  /// Configuration containing the finding function instance
  ActsExamples::TrkrClusterFindingAlgorithm::Config findCfg;

  Acts::PropagatorPlainOptions m_actsPropPlainOptions;
};

#endif
