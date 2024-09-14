// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_TPCCLUSTERQA_H
#define QA_TRACKING_TPCCLUSTERQA_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <vector>

class PHCompositeNode;
class TH2;
class TH1;
class TpcClusterQA : public SubsysReco
{
 public:
  TpcClusterQA(const std::string &name = "TpcClusterQA");

  ~TpcClusterQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

 private:
  void createHistos();

  std::vector<float> m_clusgz;
  std::vector<int> m_cluslayer;
  std::vector<int> m_clusphisize;
  std::vector<int> m_cluszsize;
  std::vector<int> m_region;

  std::string getHistoPrefix() const;
  std::set<int> m_layers;
  std::multimap<int, int> m_layerRegionMap;

  int m_event = 0;
  int m_totalClusters = 0;
  int m_clustersPerSector[24] = {0};

  TH2 *h_totalclusters = nullptr;
  TH2 *h_clusterssector = nullptr;
  TH2 *h_hitpositions = nullptr;
  TH1 *h_hitzpositions_side0 = nullptr;
  TH1 *h_hitzpositions_side1 = nullptr;

  TH1 *h_phisize_side0[3] = {nullptr};
  TH1 *h_phisize_side1[3] = {nullptr};
  TH1 *h_zsize[3] = {nullptr};
  TH1 *h_rphierror[3] = {nullptr};
  TH1 *h_zerror[3] = {nullptr};
  TH1 *h_clusedge[3] = {nullptr};
  TH1 *h_clusoverlap[3] = {nullptr};
  TH1 *h_clusxposition_side0[3] = {nullptr};
  TH1 *h_clusxposition_side1[3] = {nullptr};
  TH1 *h_clusyposition_side0[3] = {nullptr};
  TH1 *h_clusyposition_side1[3] = {nullptr};
  TH1 *h_cluszposition_side0[3] = {nullptr};
  TH1 *h_cluszposition_side1[3] = {nullptr};
};

#endif  // QA_TRACKING_TPCCLUSTERQA_H
