// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHSILICONCOSMICSEEDING_H
#define PHSILICONCOSMICSEEDING_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <map>
#include <set>
#include <string>
#include <vector>

class PHCompositeNode;
class TrkrClusterContainer;
class ActsGeometry;
class TrackSeedContainer;
class PHSiliconCosmicSeeding : public SubsysReco
{
 public:
  using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;
  struct seed
  {
    float xyslope = std::numeric_limits<float>::quiet_NaN();
    float xyintercept = std::numeric_limits<float>::quiet_NaN();
    float rzslope = std::numeric_limits<float>::quiet_NaN();
    float rzintercept = std::numeric_limits<float>::quiet_NaN();
    std::set<TrkrDefs::cluskey> ckeys;
  };
  using SeedVector = std::vector<seed>;
  PHSiliconCosmicSeeding(const std::string &name = "PHSiliconCosmicSeeding");

  ~PHSiliconCosmicSeeding() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  SeedVector makeDoublets(PositionMap &clusterPositions);
  SeedVector addClustersOnLine(SeedVector &doublets, PositionMap &clusterPositions);
  SeedVector combineSeeds(SeedVector &doublets);
  void pruneSeeds(SeedVector &doublets, PositionMap &clusterPositions);
  TrkrClusterContainer *m_clusterContainer = nullptr;
  ActsGeometry *m_tGeometry = nullptr;
  TrackSeedContainer *m_seedContainer = nullptr;
  float m_maxDoubletDistance = 6.;  // cm
  std::string m_trackMapName = "SiliconTrackSeedContainer";
};

#endif  // PHSILICONCOSMICSEEDING_H
