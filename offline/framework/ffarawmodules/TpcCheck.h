// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_TPCCHECK_H
#define FFARAWMODULES_TPCCHECK_H

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class TpcCheck : public SubsysReco
{
 public:
  TpcCheck(const std::string &name = "TpcCheck");

  ~TpcCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

//  int ResetEvent(PHCompositeNode *topNode) override;

  void MyEvtNode(const std::string &name) {m_EvtNodeName = name;}

  void SetBcoRange(const unsigned int i) {bcorange = i;}
 private:
  unsigned int bcorange {0};
  std::string m_EvtNodeName = "TPCRAWHIT";
  std::set<uint64_t> bclk_seen;
};

#endif  // FFARAWMODULES_EVENTCOMBINER_H
