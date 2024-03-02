// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef MVTXCLUSTERQA_H
#define MVTXCLUSTERQA_H

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

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;

  void beginRun(const int run) { m_beginRun = run; }
  void endRun(const int run) { m_endRun = run; }

 private:
  void createHistos();

  std::string getHistoPrefix() const;
  std::map<int, int> m_layerStaveMap;
  int m_event = 0;
  int m_totalClusters = 0;
  int m_beginRun = 25900;
  int m_endRun = 26200;
  int m_runbins = m_endRun - m_beginRun;
  int m_nclustersPerChip[3][20][9] = {{{0}}};
};

#endif  // MvtxClusterQA_H
