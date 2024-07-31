// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_MVTXCLUSTERQA_H
#define QA_TRACKING_MVTXCLUSTERQA_H

#include <fun4all/SubsysReco.h>

#include <cmath>
#include <map>
#include <set>
#include <string>

class PHCompositeNode;
class TH1;
class TH2;

class MvtxClusterQA : public SubsysReco
{
 public:
  MvtxClusterQA(const std::string &name = "MvtxClusterQA");

  ~MvtxClusterQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

  void writeChipInfo(bool value)
  {
    m_chipInfo = value;
  }

 private:
  void createHistos();

  std::string getHistoPrefix() const;
  std::map<int, int> m_layerStaveMap;
  int m_event = 0;
  int m_totalClusters = 0;
  int m_nclustersPerChip[3][20][9] = {{{0}}};
  bool m_chipInfo = false;

  TH1 *h_occupancy{nullptr};
  TH1 *h_clusSize{nullptr};
  TH2 *h_clusSize_nClus{nullptr};
  TH1 *h_clusPhi_incl{nullptr};
  TH1 *h_clusPhi_l0{nullptr};
  TH1 *h_clusPhi_l1{nullptr};
  TH1 *h_clusPhi_l2{nullptr};
  TH2 *h_clusZ_clusPhi_l0{nullptr};
  TH2 *h_clusZ_clusPhi_l1{nullptr};
  TH2 *h_clusZ_clusPhi_l2{nullptr};
  TH1 *h_strobe{nullptr};
  TH2 *h_clusperchip[3][20][9] = {{{nullptr}}};
};

#endif  // QA_TRACKING_MVTXCLUSTERQA_H
