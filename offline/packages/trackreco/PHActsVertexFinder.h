#ifndef TRACKRECO_PHACTSVERTEXFINDER_H
#define TRACKRECO_PHACTSVERTEXFINDER_H

#include "PHInitVertexing.h"
#include "ActsTrackingGeometry.h"

#include <trackbase/TrkrDefs.h>

#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;

namespace Acts
{
  class TrackParameters;
}

using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;


class PHActsVertexFinder: public PHInitVertexing 
{
  
 public:
  PHActsVertexFinder(const std::string &name);
  virtual ~PHActsVertexFinder() {}



 protected:
  int Setup(PHCompositeNode *topNode) override;
  int Process(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
 private:
  
  int createNodes(PHCompositeNode *topNode);
  int getNodes(PHCompositeNode *topNode);
  std::vector<const Acts::BoundTrackParameters*> getTracks();
  std::map<const unsigned int, Trajectory> *m_actsFitResults;

  int m_event;
  int m_maxVertices;
  ActsTrackingGeometry *m_tGeometry;
  
  
};

#endif // TRACKRECO_PHACTSVERTEXFINDER_H
