#ifndef TRACKRECO_PHACTSINITIALVERTEXFINDER_H
#define TRACKRECO_PHACTSINITIALVERTEXFINDER_H

#include "PHInitVertexing.h"
#include <trackbase/ActsTrackingGeometry.h>

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

using InitKeyMap = std::map<const Acts::BoundTrackParameters*, const unsigned int>;

class PHActsInitialVertexFinder: public PHInitVertexing
{
 public: 
  PHActsInitialVertexFinder(const std::string& name="PHActsInitialVertexFinder");
  virtual ~PHActsInitialVertexFinder() {}

  void setMaxVertices(const int maxVertices)
  { m_maxVertices = maxVertices;}

  void setSvtxVertexMapName(const std::string& name)
  { m_svtxVertexMapName = name; }
  
  void setSvtxTrackMapName(const std::string& name)
  { m_svtxTrackMapName = name; }

  void setInitialVertexer(const bool initial)
  { m_initial = initial; }

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
  void createDummyVertex();
  void checkTrackVertexAssociation();
  
  int m_maxVertices = 5;
  int m_event = 0;
  unsigned int m_totVertexFits = 0;
  unsigned int m_successFits = 0;

  std::string m_svtxTrackMapName = "SvtxSiliconTrackMap";
  std::string m_svtxVertexMapName = "SvtxVertexMap";
  bool m_initial = true;

  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  ActsTrackingGeometry *m_tGeometry = nullptr;

};


#endif
