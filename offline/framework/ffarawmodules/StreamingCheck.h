// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_STREAMINGCHECKCHECK_H
#define FFARAWMODULES_STREAMINGCHECKCHECK_H

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class StreamingCheck : public SubsysReco
{
 public:
  StreamingCheck(const std::string &name = "StreamingCheck");

  ~StreamingCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

//  int ResetEvent(PHCompositeNode *topNode) override;

  void SetTpcBcoRange(const unsigned int i) {tpc_bcorange = i;}
 private:
  unsigned int tpc_bcorange {0};
  std::set<uint64_t> bclk_seen;
};

#endif  // FFARAWMODULES_EVENTCOMBINER_H
