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
class TpcDistortionCorrectionContainer;
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

  void makeResidQAHistos(bool value)
  {
    m_residQA = value;
  }

 private:
  void createHistos();

  bool m_residQA = false;

  std::string m_trackMapName = "SvtxTrackMap";
  TpcDistortionCorrectionContainer *m_dccModuleEdge{nullptr}, *m_dccStatic{nullptr}, *m_dccAverage{nullptr}, *m_dccFluctuation{nullptr};
  float m_px = std::numeric_limits<float>::quiet_NaN();
  float m_py = std::numeric_limits<float>::quiet_NaN();
  float m_pt = std::numeric_limits<float>::quiet_NaN();
  int m_ntpc = std::numeric_limits<int>::quiet_NaN();
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
  TH1 *h_ntpc = nullptr;

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
  TH1 *h_clusphisize1pt_side0[3] = {nullptr};
  TH1 *h_clusphisize1pt_side1[3] = {nullptr};
  TH1 *h_clusphisizegeq1pt_side0[3] = {nullptr};
  TH1 *h_clusphisizegeq1pt_side1[3] = {nullptr};
};

#endif  // QA_TRACKING_TPCCLUSTERQA_H
