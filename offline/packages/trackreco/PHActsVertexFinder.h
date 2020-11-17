#ifndef TRACKRECO_PHACTSVERTEXFINDER_H
#define TRACKRECO_PHACTSVERTEXFINDER_H

#include "PHInitVertexing.h"
#include "ActsTrackingGeometry.h"

#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/Result.hpp>
#include <Acts/Vertexing/Vertex.hpp>

#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;

namespace Acts
{
  class TrackParameters;
}

using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;
using VertexVector = std::vector<Acts::Vertex<Acts::BoundTrackParameters>>;
using TrackPtrVector = std::vector<const Acts::BoundTrackParameters*>;
using VertexMap = std::map<unsigned int, Acts::Vertex<Acts::BoundTrackParameters>>;


class PHActsVertexFinder: public PHInitVertexing 
{
  
 public:
  PHActsVertexFinder(const std::string &name = "PHActsVertexFinder");
  virtual ~PHActsVertexFinder() {}
  void setMaxVertices(int maxVertices)
  { m_maxVertices = maxVertices; }

 protected:
  int Setup(PHCompositeNode *topNode) override;
  int Process(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  
  int createNodes(PHCompositeNode *topNode);
  int getNodes(PHCompositeNode *topNode);
  std::vector<const Acts::BoundTrackParameters*> getTracks();
  std::map<const unsigned int, Trajectory> *m_actsFitResults;
  VertexVector findVertices(TrackPtrVector& tracks);
  void fillVertexMap(VertexVector& vertices);
  
  int m_event = 0;
  int m_maxVertices = 15;
  VertexMap *m_vertexMap = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;
    
};

#endif // TRACKRECO_PHACTSVERTEXFINDER_H
