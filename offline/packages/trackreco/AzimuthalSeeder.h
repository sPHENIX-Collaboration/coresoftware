// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef AZIMUTHALSEEDER_H
#define AZIMUTHALSEEDER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>
class PHCompositeNode;

class ActsGeometry;
class TrkrClusterContainer;
class TrackSeedContainer;
class TFile;
class TH2;
class AzimuthalSeeder : public SubsysReco
{
 public:
  using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;
  struct seed
  {
    std::vector<TrkrDefs::cluskey> ckeys;
  };
  using SeedVector = std::vector<seed>;
  AzimuthalSeeder(const std::string &name = "AzimuthalSeeder");

  ~AzimuthalSeeder() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  TFile *file = nullptr;
  TH2 *h_phi = nullptr;
  TH2 *h_phi2 = nullptr;
  TH2 *h_phi3 = nullptr;

  std::string m_trackMapName = "SiliconTrackSeedContainer";
  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  TrackSeedContainer *m_seedContainer = nullptr;
};

#endif  // AZIMUTHALSEEDER_H
