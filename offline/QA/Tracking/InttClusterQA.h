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

  void beginRun(const int run) { m_beginRun = run; }
  void endRun(const int run) { m_endRun = run; }

 private:
  void createHistos();

  std::string getHistoPrefix() const;
  std::map<int, int> m_layerLadderMap;
  int m_event = 0;
  int m_totalClusters = 0;
  int m_beginRun = 25900;
  int m_endRun = 26200;
  int m_runbins = m_endRun - m_beginRun;
  int m_nclustersPerSensor[4][16][4] = {{{0}}};
};

#endif  // InttClusterQA_H
