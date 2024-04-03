// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHCOSMICSEEDCOMBINER_H
#define PHCOSMICSEEDCOMBINER_H

#include <fun4all/SubsysReco.h>

#include <string>

class ActsGeometry;
class TrkrClusterContainer;
class TrackSeedContainer;
class PHCompositeNode;
class TrackSeed;

class PHCosmicSeedCombiner : public SubsysReco
{
 public:
  PHCosmicSeedCombiner(const std::string &name = "PHCosmicSeedCombiner");

  ~PHCosmicSeedCombiner() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void dphiCut(float dphi) { m_dphiCut = dphi; }
  void detaCut(float deta) { m_detaCut = deta; }

 private:
  int getNodes(PHCompositeNode *topNode);
  void addKeys(TrackSeed *seedToAddTo, TrackSeed *seedToAdd);

  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  TrackSeedContainer *m_seedMap = nullptr;
  TrackSeedContainer *m_tpcSeeds = nullptr;
  TrackSeedContainer *m_siliconSeeds = nullptr;

  float m_dphiCut = 0.04;
  float m_detaCut = 0.02;
};

#endif  // PHCOSMICSEEDCOMBINER_H
