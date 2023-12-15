// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHCOSMICSEEDER_H
#define PHCOSMICSEEDER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrDefs.h>

#include <map>
#include <string>
class PHCompositeNode;
class ActsGeometry;
class TrkrClusterContainer;
class TrackSeedContainer;
class TFile;
class TNtuple;
class PHCosmicSeeder : public SubsysReco
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
  PHCosmicSeeder(const std::string &name = "PHCosmicSeeder");

  ~PHCosmicSeeder() override;
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void xyTolerance(float tol) { m_xyTolerance = tol; }
  void seedAnalysis() { m_analysis = true; }

 private:
  int getNodes(PHCompositeNode *topNode);
  int createNodes(PHCompositeNode *topNode);
  SeedVector makeSeeds(PositionMap &clusterPositions);
  SeedVector combineSeeds(SeedVector &initialSeeds, PositionMap &clusterPositions);
  SeedVector findIntersections(SeedVector &initialSeeds);
  SeedVector chainSeeds(SeedVector &initialSeeds, PositionMap &clusterPositions);
  void recalculateSeedLineParameters(seed &seed, PositionMap &clusters, bool isXY);

  float m_xyTolerance = 2.;  //! cm
  float m_rzTolerance = 2.;  //! cm
  std::string m_trackMapName = "TpcTrackSeedContainer";
  ActsGeometry *m_tGeometry = nullptr;
  TrkrClusterContainer *m_clusterContainer = nullptr;
  TrackSeedContainer *m_seedContainer = nullptr;
  TFile *m_outfile = nullptr;
  TNtuple *m_tup = nullptr;
  bool m_analysis = false;
  float m_event = 0;
};

#endif  // PHCOSMICSEEDER_H
