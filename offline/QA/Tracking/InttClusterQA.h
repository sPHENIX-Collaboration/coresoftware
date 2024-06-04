// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef INTTCLUSTERQA_H
#define INTTCLUSTERQA_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>

class PHCompositeNode;

class InttClusterQA : public SubsysReco
{
 public:
  InttClusterQA(const std::string &name = "InttClusterQA");

  ~InttClusterQA() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;

  void writeSensorInfo(bool value)
  {
    m_sensorInfo = value;
  }

 private:
  void createHistos();

  std::string getHistoPrefix() const;
  std::map<int, int> m_layerLadderMap;
  int m_event = 0;
  int m_totalClusters = 0;
  int m_nclustersPerSensor[4][16][4] = {{{0}}};
  bool m_sensorInfo = false;
};

#endif  // InttClusterQA_H
