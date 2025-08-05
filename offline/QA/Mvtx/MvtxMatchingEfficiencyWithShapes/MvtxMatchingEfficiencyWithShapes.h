// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_MVTXMATHINGEFFICIENCYWITHSHAPES_H_
#define QA_TRACKING_MVTXMATHINGEFFICIENCYWITHSHAPES_H_

#include <fun4all/SubsysReco.h>

#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <trackbase/TrkrDefs.h>

class ActsGeometry;
class TrkrHitSetContainer;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class PHG4TpcCylinderGeomContainer;
class TrkrHitSetContainer;
class TrkrClusterHitAssoc;


class TH1;
class TH1D;
class TH2D;
class TH1I;
class TH2;
class TCanvas;
class TrackSeed;
class TTree;

class PHCompositeNode;

class MvtxMatchingEfficiencyWithShapes : public SubsysReco
{
 public:
 MvtxMatchingEfficiencyWithShapes(const std::string &name = "MvtxMatchingEfficiencyWithShapes",const std::string outputfilename = "out.root");

  ~MvtxMatchingEfficiencyWithShapes() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

  float calc_dedx(TrackSeed *tpcseed);
  std::set<TrkrDefs::cluskey> findDuplicates(std::vector<TrkrDefs::cluskey> vec);

  private:

  std::string m_outputFileName;

  float pt, eta, phi, frac_p_z, dEdx;
  int layers, states,nTPC;

  int ievent = 0;
  TH1I *h_status = nullptr;
  TH1I *h_INTT_time_delta = nullptr;

  TrkrClusterContainer *_cluster_map = nullptr;
  PHG4TpcCylinderGeomContainer *_geom_container{nullptr};
  ActsGeometry *m_tGeometry = nullptr;

  //! hits
  TrkrHitSetContainer* m_hitsetcontainer = nullptr;
  //! cluster to hit association
  TrkrClusterHitAssoc* m_cluster_hit_map = nullptr;

static constexpr int nLayers = 3;
static constexpr int nClustersPerLayer = 2;

int nhits_arr[nLayers][nClustersPerLayer] = {};
int x_arr[nLayers][nClustersPerLayer] = {};
int y_arr[nLayers][nClustersPerLayer] = {};
unsigned long long key_arr[nLayers][nClustersPerLayer] = {};

  TTree *tree = nullptr;
};

#endif  // QA_TRACKING_MVTXMATHINGEFFICIENCYWITHSHAPES_H_
