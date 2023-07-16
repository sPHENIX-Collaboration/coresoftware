// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_EVTCHECK_H
#define FFARAWMODULES_EVTCHECK_H

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class EvtCheck : public SubsysReco
{
 public:
  EvtCheck(const std::string &name = "EvtCheck");

  ~EvtCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

//  int ResetEvent(PHCompositeNode *topNode) override;

  void MyEvtNode(const std::string &name) {m_EvtNodeName = name;}

 private:
  std::string m_EvtNodeName = "EVT";
};

#endif  // FFARAWMODULES_EVENTCOMBINER_H
