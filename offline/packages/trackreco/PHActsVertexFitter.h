#ifndef TRACKRECO_PHACTSVERTEXFITTER_H
#define TRACKRECO_PHACTSVERTEXFITTER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsTrackingGeometry.h>

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

using BoundTrackParamVec = std::vector<const Acts::BoundTrackParameters *>;
using VertexTrackMap = std::map<const unsigned int,
                                BoundTrackParamVec>;

using ActsVertex = const Acts::Vertex<Acts::BoundTrackParameters>;

/**
 * This class runs the Acts vertex fitter on the final tracks. It is
 * required that the tracks already have an identified vertexId associated
 * to them, i.e. that vertex finding has already been performed
 */
class PHActsVertexFitter : public SubsysReco
{
 public:
  PHActsVertexFitter(const std::string &name = "PHActsVertexFitter");
  ~PHActsVertexFitter() override {}
  int process_event(PHCompositeNode *topNode) override;
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void updateSvtxVertexMap(bool updateSvtxVertexMap)
  {
    m_updateSvtxVertexMap = updateSvtxVertexMap;
  }

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);

  /// Get the tracks with their associated vertex Ids
  VertexTrackMap getTracks();

  /// Turn the SvtxTrack object into an Acts::TrackParameters object
  const Acts::BoundTrackParameters *makeTrackParam(const SvtxTrack *track) const;

  /// Run the Acts vertex fitter
  ActsVertex fitVertex(BoundTrackParamVec tracks,
                       Acts::Logging::Level logLevel) const;

  /// Runs Acts vertex fitter
  void fitVertices(std::vector<const Acts::BoundTrackParameters *> tracks);

  /// Update SvtxVertex or create new SvtxVertexMap
  void createActsSvtxVertex(const unsigned int,
                            ActsVertex vertex);
  void updateSvtxVertex(const unsigned int,
                        ActsVertex vertex);

  int m_event = 0;

  std::map<const unsigned int, Trajectory> *m_actsFitResults;
  ActsTrackingGeometry *m_tGeometry;
  SvtxTrackMap *m_trackMap = nullptr;
  SvtxVertexMap *m_vertexMap = nullptr;
  SvtxVertexMap *m_actsVertexMap = nullptr;

  /// Option to update the default SvtxVertexMap. A new SvtxVertexMap
  /// called SvtxVertexMapActs is created by default in the module
  bool m_updateSvtxVertexMap = false;
};

#endif  //TRACKRECO_PHACTSVERTEXFITTER_H
