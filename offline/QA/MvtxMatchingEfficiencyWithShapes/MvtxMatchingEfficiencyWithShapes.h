// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_MVTXMATHINGEFFICIENCYWITHSHAPES_H_
#define QA_TRACKING_MVTXMATHINGEFFICIENCYWITHSHAPES_H_

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>

#include <array>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <vector>

class ActsGeometry;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class PHG4TpcGeomContainer;
class TrkrClusterHitAssoc;
class TH1;
class TrackSeed;
class TTree;
class PHCompositeNode;

class MvtxMatchingEfficiencyWithShapes : public SubsysReco
{
 public:
  MvtxMatchingEfficiencyWithShapes(const std::string &name = "MvtxMatchingEfficiencyWithShapes", const std::string &outputfilename = "out.root");

  ~MvtxMatchingEfficiencyWithShapes() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

  float calc_dedx(TrackSeed *tpcseed);
  static std::set<TrkrDefs::cluskey> findDuplicates(std::vector<TrkrDefs::cluskey> vec);

 private:
  std::string m_outputFileName;

  float pt{std::numeric_limits<float>::quiet_NaN()};
  float eta{std::numeric_limits<float>::quiet_NaN()};
  float phi{std::numeric_limits<float>::quiet_NaN()};
  float frac_p_z{std::numeric_limits<float>::quiet_NaN()};
  float dEdx{std::numeric_limits<float>::quiet_NaN()};
  int layers{std::numeric_limits<int>::min()};
  int states{std::numeric_limits<int>::min()};
  int nTPC{std::numeric_limits<int>::min()};

  int ievent {0};
  TH1 *h_status {nullptr};
  TH1 *h_INTT_time_delta {nullptr};

  TrkrClusterContainer *_cluster_map {nullptr};
  PHG4TpcGeomContainer *_geom_container{nullptr};
  ActsGeometry *m_tGeometry {nullptr};

  //! hits
  TrkrHitSetContainer *m_hitsetcontainer {nullptr};
  //! cluster to hit association
  TrkrClusterHitAssoc *m_cluster_hit_map {nullptr};

  static constexpr int nLayers {3};
  static constexpr int nClustersPerLayer {2};
  std::array<std::array<int, nClustersPerLayer>, nLayers> nhits_arr{};
  std::array<std::array<int, nClustersPerLayer>, nLayers> x_arr{};
  std::array<std::array<int, nClustersPerLayer>, nLayers> y_arr{};
  std::array<std::array<uint64_t, nClustersPerLayer>, nLayers> key_arr{};

  TTree *tree {nullptr};
};

#endif  // QA_TRACKING_MVTXMATHINGEFFICIENCYWITHSHAPES_H_
