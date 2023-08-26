// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHCOSMICSEEDCOMBINER_H
#define PHCOSMICSEEDCOMBINER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <TH1.h>
#include <TFile.h>

class ActsGeometry;
class TrkrClusterContainer;
class TrackSeedContainer;
class PHCompositeNode;

class PHCosmicSeedCombiner : public SubsysReco
{
 public:

  PHCosmicSeedCombiner(const std::string &name = "PHCosmicSeedCombiner");

  ~PHCosmicSeedCombiner() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);

  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  TrackSeedContainer *m_seedMap = nullptr;
  TrackSeedContainer *m_tpcSeeds = nullptr;
  TrackSeedContainer *m_siliconSeeds = nullptr;
  TrackSeedContainer *m_cosmicContainer = nullptr;
  
  TFile *file = nullptr;
  TH1 *h_deta = nullptr;
  TH1 *h_dphi = nullptr;
};

#endif // PHCOSMICSEEDCOMBINER_H
