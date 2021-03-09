#ifndef TRACKRECO_PHACTSTOSVTXTRACKS_H
#define TRACKRECO_PHACTSTOSVTXTRACKS_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>

/// Acts includes to create all necessary definitions
#include <Acts/Utilities/BinnedArray.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <Acts/EventData/TrackParameters.hpp>

#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>
#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

#include <trackbase/ActsTrackingGeometry.h>

#include <boost/bimap.hpp>

#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class MakeActsGeometry;

using SourceLink = ActsExamples::TrkrClusterSourceLink;
using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;

typedef boost::bimap<TrkrDefs::cluskey, unsigned int> CluskeyBimap;


class PHActsToSvtxTracks : public SubsysReco
{
 public:
  /// Default constructor and destructor
  PHActsToSvtxTracks(const std::string &name = "PHActsToSvtxTracks");
  virtual ~PHActsToSvtxTracks() {}

  /// Inherited SubsysReco functions
  int End(PHCompositeNode *topNode);
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  
  void setSvtxTrackMapName(std::string &name)
  { m_svtxMapName = name;}

 private:

  /// Create track seed node if it doesn't exist yet
  int createNodes(PHCompositeNode *topNode);

  /// Get nodes off node tree needed to execute module
  int getNodes(PHCompositeNode *topNode);

  void createSvtxTrack(const unsigned int trackKey, Trajectory traj);

  SvtxTrackMap *m_svtxTrackMap = nullptr;
  SvtxVertexMap *m_svtxVertexMap = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;
  std::map<const unsigned int, Trajectory> *m_actsFitResults = nullptr;
  CluskeyBimap *m_hitIdClusKey = nullptr;
  std::string m_svtxMapName = "SvtxTrackMap";

};

#endif
