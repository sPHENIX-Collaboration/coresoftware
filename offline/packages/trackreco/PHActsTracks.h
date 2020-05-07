#ifndef TRACKRECO_PHACTSTRACKS_H
#define TRACKRECO_PHACTSTRACKS_H

#include "PHActsSourceLinks.h" 

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/TrackParameters.hpp>

#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/EventData/TrkrClusterSourceLink.hpp>

#include "ActsTrack.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class MakeActsGeometry;

using SourceLink = FW::Data::TrkrClusterSourceLink;


/**
 * This class is responsible for taking SvtxTracks and converting them to track
 * seeds that Acts can take in to the track fitter. It collects SvtxTracks, 
 * converts them to Acts tracks, and then finds the corresponding
 * TrkrClusterSourceLinks to that SvtxTrack. The output is a node on the node
 * tree that is a map of Acts track seeds and corresponding source links
 */
class PHActsTracks : public SubsysReco
{
 public:
  /// Default constructor and destructor
  PHActsTracks(const std::string &name = "PHActsTracks");
  virtual ~PHActsTracks() {}

  /// Inherited SubsysReco functions
  int End(PHCompositeNode *topNode);
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);

 private:
  /** 
   * Member functions
   */

  /// Create track seed node if it doesn't exist yet
  void createNodes(PHCompositeNode *topNode);

  /// Get nodes off node tree needed to execute module
  int getNodes(PHCompositeNode *topNode);

  Acts::BoundSymMatrix getActsCovMatrix(const SvtxTrack *track);

  /**
   * Member variables
   */

  TFile *outfile;
  TH1F *h_x, *h_y, *h_z;
  /// A vector to hold the source links corresponding to a particular SvtxTrack
  std::vector<SourceLink> m_trackSourceLinks;

  /// A map of an Acts track seed and Acts-sPHENIX source links corresponding
  /// to that track seed
  std::vector<ActsTrack> *m_actsProtoTracks;

  /// Trackmap that contains SvtxTracks
  SvtxTrackMap *m_trackMap;

  /// Map between cluster key and arbitrary hit id created in PHActsSourceLinks
  std::map<TrkrDefs::cluskey, unsigned int> *m_hitIdClusKey;

  /// Map of hitid:SourceLinks created in PHActsSourceLinks
  std::map<unsigned int, SourceLink> *m_sourceLinks;

  ActsTrackingGeometry *m_tGeometry;
};

#endif
