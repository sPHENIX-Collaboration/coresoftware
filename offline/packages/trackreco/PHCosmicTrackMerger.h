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

  void zero_field() { m_zeroField = true; }
  void dca_xycut(const float cut) { m_dcaxycut = cut; }
  void dca_rzcut(const float cut) { m_dcarzcut = cut; }
 private:
  void addKeys(TrackSeed *toAddTo, TrackSeed *toAdd);
  void removeOutliers(TrackSeed *seed);

  ActsGeometry *m_geometry = nullptr;
  KeyPosMap getGlobalPositions(TrackSeed *seed);
  TrkrClusterContainer *m_clusterMap = nullptr;
  TrackSeedContainer *m_seeds = nullptr;
  TrackSeedContainer *m_tpcSeeds = nullptr;
  TrackSeedContainer *m_siliconSeeds = nullptr;

  float m_dcaxycut = 0.5; // cm
  float m_dcarzcut = 2.; // cm
  bool m_zeroField = false;
};

#endif  // PHCOSMICTRACKMERGER_H
