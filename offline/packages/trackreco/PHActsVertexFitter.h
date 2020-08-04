#ifndef TRACKRECO_PHACTSVERTEXFITTER_H
#define TRACKRECO_PHACTSVERTEXFITTER_H

#include <fun4all/SubsysReco.h>
#include "PHActsSourceLinks.h"

#include <ACTFW/EventData/TrkrClusterMultiTrajectory.hpp>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;

namespace Acts
{
class TrackParameters;
}

using Trajectory = FW::TrkrClusterMultiTrajectory;

class PHActsVertexFitter : public SubsysReco
{
 public:
  PHActsVertexFitter();
  virtual ~PHActsVertexFitter(){}
  int process_event(PHCompositeNode *topNode);
  int Init(PHCompositeNode *topNode);
  int End (PHCompositeNode *topNode);


 private:

  int getNodes(PHCompositeNode *topNode);
std::vector<Acts::BoundParameters> getTracks();  
  std::map<const unsigned int, Trajectory> *m_actsFitResults;

  ActsTrackingGeometry *m_tGeometry;
};

#endif //TRACKRECO_PHACTSVERTEXFITTER_H 
