// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHCOSMICSEEDER_H
#define PHCOSMICSEEDER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>

#include <map>
#include <string>
class PHCompositeNode;
class ActsGeometry;
class TrkrClusterContainer;
class TrackSeedContainer;

class PHCosmicSeeder : public SubsysReco
{
 public:
  using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;
  struct seed
  {
    float xyslope;
    float xyintercept;
    float rzslope;
    float rzintercept;
    std::vector<TrkrDefs::cluskey> ckeys;
  };
  using SeedVec = std::vector<seed>;
  PHCosmicSeeder(const std::string &name = "PHCosmicSeeder");

  ~PHCosmicSeeder() override;
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  SeedVec makeSeeds(PositionMap &clusterPositions);

  

  float m_xyTolerance = 2.; //! cm
  float m_rzTolerance = 2.; //! cm
  std::string m_trackMapName = "TpcTrackSeedContainer";
  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  TrackSeedContainer *m_seedContainer = nullptr;
};

#endif // PHCOSMICSEEDER_H
