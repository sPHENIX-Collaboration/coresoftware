
#ifndef TRACKSEEDTRACKMAPCONVERTER_H
#define TRACKSEEDTRACKMAPCONVERTER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <memory>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack_v4;
class TrackSeed;
class TrackSeedContainer;
class ActsGeometry;
class TrkrClusterContainer;

class TrackSeedTrackMapConverter : public SubsysReco
{
 public:

  TrackSeedTrackMapConverter(const std::string &name = "TrackSeedTrackMapConverter");

  virtual ~TrackSeedTrackMapConverter();

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void setTrackMapName(const std::string& name) { m_trackMapName = name; }
  void setTrackSeedName(const std::string& name) { m_trackSeedName = name; }

 private:

  int getNodes(PHCompositeNode *topNode);

  void addKeys(std::unique_ptr<SvtxTrack_v4>& track, TrackSeed *seed);
  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_trackSeedName = "TpcTrackSeedContainer";

  SvtxTrackMap *m_trackMap = nullptr;
  TrackSeedContainer *m_seedContainer = nullptr;
  TrackSeedContainer *m_tpcContainer = nullptr;
  TrackSeedContainer *m_siContainer = nullptr;

  TrkrClusterContainer *m_clusters = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
};

#endif // TRACKSEEDTRACKMAPCONVERTER_H
