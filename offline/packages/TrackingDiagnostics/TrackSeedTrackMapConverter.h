
#ifndef TRACKSEEDTRACKMAPCONVERTER_H
#define TRACKSEEDTRACKMAPCONVERTER_H

#include <fun4all/SubsysReco.h>

#include <limits>
#include <memory>
#include <string>

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

  virtual ~TrackSeedTrackMapConverter() = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void setFieldMap(const std::string &name) { m_fieldMap = name; }
  void setTrackMapName(const std::string &name) { m_trackMapName = name; }
  void setTrackSeedName(const std::string &name) { m_trackSeedName = name; }
  void cosmics() { m_cosmics = true; }

 private:
  int getNodes(PHCompositeNode *topNode);

  void addKeys(std::unique_ptr<SvtxTrack_v4> &track, TrackSeed *seed);
  void addKeys(TrackSeed *seedToAddTo, TrackSeed *seedToAdd);
  std::pair<int, float> getCosmicCharge(TrackSeed *seed, float vertexradius) const;

  ActsGeometry *m_tGeometry{nullptr};
  SvtxTrackMap *m_trackMap{nullptr};
  TrackSeedContainer *m_seedContainer{nullptr};
  TrackSeedContainer *m_tpcContainer{nullptr};
  TrackSeedContainer *m_siContainer{nullptr};
  TrkrClusterContainer *m_clusters{nullptr};

  double fieldstrength{std::numeric_limits<double>::quiet_NaN()};

  bool m_cosmics{false};
  bool m_ConstField{false};

  std::string m_fieldMap;
  std::string m_trackMapName{"SvtxTrackMap"};
  std::string m_trackSeedName{"TpcTrackSeedContainer"};
};

#endif  // TRACKSEEDTRACKMAPCONVERTER_H
