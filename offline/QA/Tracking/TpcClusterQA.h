// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_TPCCLUSTERQA_H
#define QA_TRACKING_TPCCLUSTERQA_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>

class PHCompositeNode;

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
    

  std::string getHistoPrefix() const;
  std::set<int> m_layers;
  std::multimap<int, int> m_layerRegionMap;
   
  int m_event = 0;
  int m_totalClusters = 0;
  int m_clustersPerSector[24] = {0};
};

#endif  // QA_TRACKING_TPCCLUSTERQA_H
