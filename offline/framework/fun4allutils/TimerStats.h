// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLUTILS_TIMERSTATS_H
#define FUN4ALLUTILS_TIMERSTATS_H

#include <fun4all/SubsysReco.h>

#include <string>  // for string

class PHCompositeNode;
class CDBTTree;

class TimerStats : public SubsysReco
{
 public:
  TimerStats(const std::string &name = "TimerStats");
  ~TimerStats() override {}
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void OutFileName(const std::string &name) { outfilename = name; }

 private:
  CDBTTree *cdbttree = nullptr;
  int iev = 0;
  std::string outfilename = "timerstats.root";
};

#endif
