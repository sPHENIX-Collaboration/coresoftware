// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_MVTXCLUSTERQA_H
#define QA_TRACKING_MVTXCLUSTERQA_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>

class PHCompositeNode;

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
};

#endif  // QA_TRACKING_MVTXCLUSTERQA_H
