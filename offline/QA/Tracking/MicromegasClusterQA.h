// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_TRACKING_MICROMEGASCLUSTERQA_H_
#define QA_TRACKING_MICROMEGASCLUSTERQA_H_

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>

class PHCompositeNode;

class MicromegasClusterQA : public SubsysReco
{
 public:
  MicromegasClusterQA(const std::string &name = "MicromegasClusterQA");

  ~MicromegasClusterQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

 private:
  void createHistos();

  std::string getHistoPrefix() const;
  std::map<int, int> m_layerTileMap;
  int m_event = 0;
  int m_totalClusters = 0;
  int m_nclustersPerTile[2][8] = {{0}};
};

#endif  // MicromegasClusterQA_H
