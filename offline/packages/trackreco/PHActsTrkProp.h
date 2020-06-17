#ifndef TRACKRECO_PHACTSTRKPROP_H
#define TRACKRECO_PHACTSTRKPROP_H

#include "PHTrackPropagating.h"
#include "PHActsSourceLinks.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/TrackFinder/CKFSourceLinkSelector.hpp>
#include <Acts/Propagator/detail/SteppingLogger.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>

#include <ACTFW/EventData/TrkrClusterSourceLink.hpp>
#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/TrackFinding/TrkrClusterFindingAlgorithm.hpp>

#include <memory>
#include <string>
#include <map>

class ActsTrack;
class MakeActsGeometry;
class SvtxTrack;
class SvtxTrackMap;

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

 private:
  /// Event counter
  int m_event;

  /// Get all the nodes
  int getNodes(PHCompositeNode *topNode);

  /// Create new nodes
  void createNodes(PHCompositeNode *topNode);

  ActsTrackingGeometry *m_tGeometry;

  /// Track map with Svtx objects
  SvtxTrackMap *m_trackMap;

  /// Acts proto tracks to be put on the node tree by this module
  std::vector<ActsTrack> *m_actsProtoTracks;

  /// Map of cluster keys to hit ids, for debugging statements
  std::map<TrkrDefs::cluskey, unsigned int> *m_hitIdClusKey;

  /// Acts source links created by PHActsSourceLinks
  /// SourceLink is defined as TrkrClusterSourceLink elsewhere
  std::map<unsigned int, SourceLink> *m_sourceLinks;

  /// SourceLinkSelector to help the CKF identify which source links 
  /// may belong to a track seed based on geometry considerations
  SourceLinkSelectorConfig m_sourceLinkSelectorConfig;
 
  /// Configuration containing the finding function instance
  FW::TrkrClusterFindingAlgorithm::Config findCfg;

};

#endif
