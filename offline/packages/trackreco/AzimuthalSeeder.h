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
    std::vector<Acts::Vector3> globpos;
  };
  using SeedVector = std::vector<seed>;
  AzimuthalSeeder(const std::string &name = "AzimuthalSeeder");

  ~AzimuthalSeeder() override = default;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void histos() { m_outfile = true; }
  void residualLimit(const float limit) { m_outlierLimit = limit; }

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  std::unique_ptr<TFile> file = nullptr;
  TH2 *h_phi = nullptr;
  TH2 *h_phi2 = nullptr;
  TH2 *h_phi3 = nullptr;
  bool m_outfile = false;
  std::string m_trackMapName = "SiliconTrackSeedContainer";
  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  TrackSeedContainer *m_seedContainer = nullptr;

  float m_outlierLimit = 0.1;
};

#endif  // AZIMUTHALSEEDER_H
