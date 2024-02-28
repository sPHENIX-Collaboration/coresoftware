// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TPCCLUSTERQA_H
#define TPCCLUSTERQA_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>

class PHCompositeNode;

class TpcClusterQA : public SubsysReco
{
 public:
  TpcClusterQA(const std::string &name = "TpcClusterQA");

  ~TpcClusterQA() override;

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
  std::set<int> m_layers;
  std::multimap<int, int> m_layerRegionMap;

  int m_event = 0;
  int m_totalClusters = 0;
  int m_beginRun = 25900;
  int m_endRun = 26200;
  int m_runbins = m_endRun - m_beginRun;
  int m_clustersPerSector[24] = {0};
};

#endif  // TPCCLUSTERQA_H
