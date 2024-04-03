// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_INTTCHECK_H
#define FFARAWMODULES_INTTCHECK_H

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class InttCheck : public SubsysReco
{
 public:
  InttCheck(const std::string &name = "InttCheck");

  ~InttCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  //  int ResetEvent(PHCompositeNode *topNode) override;

  void MyEvtNode(const std::string &name) { m_EvtNodeName = name; }

 private:
  std::string m_EvtNodeName = "INTTRAWHIT";
  std::set<uint64_t> bclk_seen;
};

#endif  // FFARAWMODULES_EVENTCOMBINER_H
