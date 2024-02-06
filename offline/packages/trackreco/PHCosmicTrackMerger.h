// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHCOSMICTRACKMERGER_H
#define PHCOSMICTRACKMERGER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsGeometry.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <string>
#include <vector>

class PHCompositeNode;
class TrackSeedContainer;
class TrkrClusterContainer;
class TrackSeed;
class ActsGeometry;

class PHCosmicTrackMerger : public SubsysReco
{
  using KeyPosMap = std::pair<std::vector<TrkrDefs::cluskey>, std::vector<Acts::Vector3>>;

 public:
  PHCosmicTrackMerger(const std::string &name = "PHCosmicTrackMerger");

  ~PHCosmicTrackMerger() override;

  int Init(PHCompositeNode *) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *) override;

 private:
  void addKeys(TrackSeed *toAddTo, TrackSeed *toAdd);
  void removeOutliers(TrackSeed *seed);

  ActsGeometry *m_geometry = nullptr;
  KeyPosMap getGlobalPositions(TrackSeed *seed);
  TrkrClusterContainer *m_clusterMap = nullptr;
  TrackSeedContainer *m_seeds = nullptr;
  TrackSeedContainer *m_tpcSeeds = nullptr;
  TrackSeedContainer *m_siliconSeeds = nullptr;
};

#endif  // PHCOSMICTRACKMERGER_H
