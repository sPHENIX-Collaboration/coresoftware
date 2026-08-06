#pragma once

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class Tpc_PolyTrack;
class Tpc_PolyTrackContainer;
class TpcCrossingDecisionContainer;
class ActsGeometry;
class TrackSeedContainer;

class TpcPolyTrackSeedConverter : public SubsysReco
{
 public:
  explicit TpcPolyTrackSeedConverter(const std::string& name = "TpcPolyTrackSeedConverter");
  ~TpcPolyTrackSeedConverter() override = default;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setInputNodeName(const std::string& name) { m_inputNodeName = name; }
  void setOutputNodeName(const std::string& name) { m_outputNodeName = name; }
  void setCrossingDecisionNodeName(const std::string& name) { m_crossingDecisionNodeName = name; }
  void setCrossingPeriodNs(double value) { m_crossingPeriodNs = value; }

 private:
  int getNodes(PHCompositeNode*);
  int createNodes(PHCompositeNode*);
  bool isValidTrack(const Tpc_PolyTrack*) const;
  bool publishSeed(const Tpc_PolyTrack*) const;

  std::string m_inputNodeName;
  std::string m_outputNodeName;
  std::string m_crossingDecisionNodeName {"TPC_CROSSING_DECISIONS"};
  Tpc_PolyTrackContainer* m_polyTracks {nullptr};
  TpcCrossingDecisionContainer* m_crossingDecisions {nullptr};
  ActsGeometry* m_geometry {nullptr};
  TrackSeedContainer* m_trackSeeds {nullptr};
  double m_crossingPeriodNs {106.56};
};
