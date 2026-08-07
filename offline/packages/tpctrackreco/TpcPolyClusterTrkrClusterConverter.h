#pragma once

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>

#include <array>
#include <map>
#include <memory>
#include <string>

class ActsGeometry;
class PHG4TpcGeomContainer;
class PHCompositeNode;
class TpcClusterMover;
class TpcCrossingDecisionContainer;
class Tpc_PolyCluster;
class Tpc_PolyClusterContainer;
class Tpc_PolyTrack;
class Tpc_PolyTrackContainer;
class TrkrClusterContainer;

class TpcPolyClusterTrkrClusterConverter : public SubsysReco
{
 public:
  explicit TpcPolyClusterTrkrClusterConverter(const std::string& name = "TpcPolyClusterTrkrClusterConverter");
  ~TpcPolyClusterTrkrClusterConverter() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setPolyClusterNodeName(const std::string& name) { m_polyClusterNodeName = name; }
  void setPolyTrackNodeName(const std::string& name) { m_polyTrackNodeName = name; }
  void setOutputNodeName(const std::string& name) { m_outputNodeName = name; }
  void setCrossingDecisionNodeName(const std::string& name) { m_crossingDecisionNodeName = name; }
  void setCrossingPeriodNs(double value) { m_crossingPeriodNs = value; }
  void setMagneticFieldTesla(double value) { m_magneticFieldTesla = value; }

 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  void clearOutputTpcClusters();
  void buildTrackMap();
  bool isAcceptedTrack(const Tpc_PolyTrack* track) const;
  bool initializeClusterMover(PHCompositeNode* topNode);
  TrkrDefs::cluskey getClusterKey(const Tpc_PolyCluster* cluster) const;
  bool seedOutputCluster(const Tpc_PolyCluster* cluster, TrkrDefs::cluskey& cluskey);
  void buildMovedClusterMap();
  bool localFromMovedGlobal(const Tpc_PolyCluster* cluster, TrkrDefs::cluskey cluskey, const std::array<double, 3>& moved_global,
                            short crossing, float& local_x, float& local_y, unsigned short& subsurfkey,
                            unsigned long long& surface_id) const;
  bool publishCluster(const Tpc_PolyCluster* cluster, const Tpc_PolyTrack* track) const;

  std::string m_polyClusterNodeName {"TPC_POLYCLUSTERS"};
  std::string m_polyTrackNodeName {"TPC_POLYTRACKS"};
  std::string m_outputNodeName {"TRKR_CLUSTER"};
  std::string m_crossingDecisionNodeName {"TPC_CROSSING_DECISIONS"};

  Tpc_PolyClusterContainer* m_polyClusters {nullptr};
  Tpc_PolyTrackContainer* m_polyTracks {nullptr};
  TrkrClusterContainer* m_outputClusters {nullptr};
  TpcCrossingDecisionContainer* m_crossingDecisions {nullptr};
  ActsGeometry* m_geometry {nullptr};
  PHG4TpcGeomContainer* m_tpcGeomContainer {nullptr};
  std::unique_ptr<TpcClusterMover> m_clusterMover;
  std::map<unsigned int, const Tpc_PolyTrack*> m_tracksBySourceId;
  std::map<TrkrDefs::cluskey, std::array<double, 3>> m_movedGlobals;
  std::map<TrkrDefs::cluskey, unsigned short> m_seedSubSurfKeys;
  double m_crossingPeriodNs {106.56};
  double m_magneticFieldTesla {1.4};
};
