#pragma once

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class Tpc_PolyTrack;
class Tpc_PolyTrackContainer;
class Tpc_PolyTrackVertexContainer;
class PHCompositeNode;

class Tpc_PolyTrackVertexer : public SubsysReco
{
 public:
  explicit Tpc_PolyTrackVertexer(const std::string& name = "Tpc_PolyTrackVertexer");
  ~Tpc_PolyTrackVertexer() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;

  void setInputNodeName(const std::string& n) { m_inputNodeName = n; }
  void setOutputNodeName(const std::string& n) { m_outputNodeName = n; }
  void setCollisionMinClusters(unsigned int v) { m_collisionMinClusters = v; }
  void setCollisionZSeparation(double v) { m_collisionZSeparation = v; }
  void setMagneticFieldTesla(double v) { m_magneticFieldTesla = v; }

 private:
  struct TrackVertexFit
  {
    bool ok {false};
    unsigned int track_id {0};
    unsigned int source_assembled_track_id {0};
    unsigned int nclusters {0};
    double dca2d {0.0};
    double z0 {0.0};
    double pca_x {0.0};
    double pca_y {0.0};
    double pca_z {0.0};
    double pca_radius {0.0};
    double pca_phi {0.0};
  };

  struct CollisionFit
  {
    bool ok {false};
    double x {0.0};
    double y {0.0};
    double z {0.0};
    double z_rms {0.0};
    unsigned int ntracks {0};
  };

  bool getNodes(PHCompositeNode* topNode);
  bool createNodes(PHCompositeNode* topNode);
  TrackVertexFit fitTrack(const Tpc_PolyTrack* trk) const;
  CollisionFit fitCollision(const std::vector<TrackVertexFit>& tracks) const;
  std::vector<CollisionFit> fitCollisions(std::vector<TrackVertexFit> tracks) const;

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  Tpc_PolyTrackContainer* m_polyTracks {nullptr};
  Tpc_PolyTrackVertexContainer* m_vertices {nullptr};
  unsigned int m_collisionMinClusters {3};
  double m_collisionZSeparation {10.0};
  double m_magneticFieldTesla {1.4};
};
