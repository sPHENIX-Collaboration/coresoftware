#ifndef TRACKRECO_PHACTSINITIALVERTEXFINDER_H
#define TRACKRECO_PHACTSINITIALVERTEXFINDER_H

#include "PHInitVertexing.h"
#include "ActsTrackingGeometry.h"

#include <trackbase/TrkrDefs.h>

#include <Acts/Utilities/Result.hpp>
#include <Acts/Vertexing/Vertex.hpp>


#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/EventData/TrkrClusterMultiTrajectory.hpp>

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertexMap;
class SvtxVertex;

using VertexVector = std::vector<Acts::Vertex<Acts::BoundTrackParameters>>;

using TrackParamVec = std::vector<const Acts::BoundTrackParameters*>;

using InitKeyMap = std::map<const ActsExamples::TrackParameters*, const unsigned int>;

class PHActsInitialVertexFinder: public PHInitVertexing
{
 public: 
  PHActsInitialVertexFinder(const std::string& name="PHActsInitialVertexFinder");
  virtual ~PHActsInitialVertexFinder() {}

  void setMaxVertices(int maxVertices)
  { m_maxVertices = maxVertices;}

 protected:
  int Setup(PHCompositeNode *topNode) override;
  int Process(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  
  TrackParamVec getTrackPointers(InitKeyMap& keyMap);
  VertexVector findVertices(TrackParamVec& tracks);
void fillVertexMap(VertexVector& vertices, InitKeyMap& keyMap);
  
  int m_maxVertices = 10;
  int m_event = 0;
  
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;

};


#endif
