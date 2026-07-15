// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRACKRECO_PHTRUTHTRACKFITTER_H
#define TRACKRECO_PHTRUTHTRACKFITTER_H

#include <fun4all/SubsysReco.h>

#include <g4main/PHG4HitDefs.h>
#include <trackbase/TrkrDefs.h>

#include <limits>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4VtxPoint;
class SvtxTrack;
class SvtxTrackMap;
class TrackSeed;
class TrackSeedContainer;
class TrkrCluster;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitTruthAssoc;

class PHTruthTrackFitter : public SubsysReco
{
 public:
  PHTruthTrackFitter(const std::string& name = "PHTruthTrackFitter");
  ~PHTruthTrackFitter() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void setTrackMapName(const std::string& name) { m_trackMapName = name; }
  void setSvtxSeedMapName(const std::string& name) { m_svtxSeedMapName = name; }
  void setTrkrClusterContainerName(const std::string& name) { m_clusterMapName = name; }
  void setDefaultCrossing(short int crossing) { m_defaultCrossing = crossing; }
  void setPositionError(float value) { m_positionError = value; }
  void setZError(float value) { m_zError = value; }
  void setExtrapolateToClusterRadius(bool value) { m_extrapolateToClusterRadius = value; }

 private:
  int createNodes(PHCompositeNode* topNode);
  int getNodes(PHCompositeNode* topNode);

  TrackSeed* getSeed(TrackSeedContainer* container, unsigned int index) const;
  unsigned int getTruthTrackId(const TrackSeed* svtxSeed, const TrackSeed* tpcSeed, const TrackSeed* siliconSeed) const;

  std::vector<const PHG4Hit*> getTruthHits(TrkrDefs::cluskey cluskey) const;
  const PHG4Hit* getG4Hit(unsigned int trkrid, PHG4HitDefs::keytype g4hitkey) const;

  bool addStateFromCluster(SvtxTrack* track, TrkrDefs::cluskey cluskey, unsigned int truthTrackId,
                           const PHG4Particle* particle, const PHG4VtxPoint* vertex,
                           unsigned int stateIndex) const;
  float getClusterRadius(TrkrDefs::cluskey cluskey, TrkrCluster* cluster) const;
  float getPathLength(const PHG4VtxPoint* vertex, float x, float y, float z, unsigned int stateIndex) const;
  int getCharge(const PHG4Particle* particle, const TrackSeed* tpcSeed, const TrackSeed* siliconSeed) const;
  short int getCrossing(const TrackSeed* tpcSeed, const TrackSeed* siliconSeed) const;

  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_svtxSeedMapName = "SvtxTrackSeedContainer";
  std::string m_clusterMapName = "TRKR_CLUSTER";

  TrackSeedContainer* m_seedMap = nullptr;
  TrackSeedContainer* m_tpcSeeds = nullptr;
  TrackSeedContainer* m_siliconSeeds = nullptr;
  SvtxTrackMap* m_trackMap = nullptr;
  TrkrClusterContainer* m_clusterMap = nullptr;
  TrkrClusterHitAssoc* m_clusterHitMap = nullptr;
  TrkrHitTruthAssoc* m_hitTruthAssoc = nullptr;
  ActsGeometry* m_tGeometry = nullptr;
  PHG4TruthInfoContainer* m_g4TruthInfo = nullptr;

  PHG4HitContainer* m_g4HitsTpc = nullptr;
  PHG4HitContainer* m_g4HitsIntt = nullptr;
  PHG4HitContainer* m_g4HitsMvtx = nullptr;
  PHG4HitContainer* m_g4HitsMicromegas = nullptr;

  short int m_defaultCrossing = 0;
  float m_positionError = 0.005;
  float m_zError = 0.01;
  bool m_extrapolateToClusterRadius = true;

  static constexpr unsigned int m_invalidTruthTrackId = std::numeric_limits<unsigned int>::max();
};

#endif
