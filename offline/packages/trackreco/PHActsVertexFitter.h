#ifndef TRACKRECO_PHACTSVERTEXFITTER_H
#define TRACKRECO_PHACTSVERTEXFITTER_H

#include <fun4all/SubsysReco.h>
#include "ActsTrackingGeometry.h"

#include <Acts/Vertexing/Vertex.hpp>

#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertex;
class SvtxVertexMap;

namespace Acts
{
  class TrackParameters;
}

using Trajectory = ActsExamples::TrkrClusterMultiTrajectory;

using BoundTrackParamVec = std::vector<const Acts::BoundTrackParameters*>;
using VertexTrackMap = std::map<const unsigned int, 
                                BoundTrackParamVec>;

using ActsVertex = const Acts::Vertex<Acts::BoundTrackParameters>;

class PHActsVertexFitter : public SubsysReco
{
 public:
  PHActsVertexFitter(const std::string& name = "PHActsVertexFitter");
  virtual ~PHActsVertexFitter(){}
  int process_event(PHCompositeNode *topNode);
  int Init(PHCompositeNode *topNode);
  int ResetEvent(PHCompositeNode *topNode);
  int End (PHCompositeNode *topNode);


 private:
  
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  
  VertexTrackMap getTracks();  
  const Acts::BoundTrackParameters* makeTrackParam(const SvtxTrack* track) const;
  ActsVertex fitVertex(BoundTrackParamVec tracks, 
			 Acts::Logging::Level logLevel) const;
  std::map<const unsigned int, Trajectory> *m_actsFitResults;

  void fitVertices(std::vector<const Acts::BoundTrackParameters*> tracks);
  
  void createActsSvtxVertex(const unsigned int,
			    ActsVertex vertex);
  void updateSvtxVertex(const unsigned int,
			ActsVertex vertex);
 
  int m_event;
  ActsTrackingGeometry *m_tGeometry;

  SvtxTrackMap *m_trackMap;
  SvtxVertexMap *m_vertexMap;
  SvtxVertexMap *m_actsVertexMap;

  bool m_updateSvtxVertexMap = false;
};

#endif //TRACKRECO_PHACTSVERTEXFITTER_H 
