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

 private:
  void createHistos();

  std::string getHistoPrefix() const;
  std::set<int> m_layers;
  std::multimap<int, int> m_layerRegionMap;

  int m_event = 0;
  int m_totalClusters = 0;
};

#endif  // TPCCLUSTERQA_H
