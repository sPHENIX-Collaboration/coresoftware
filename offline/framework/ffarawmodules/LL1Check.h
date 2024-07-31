// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_LL1CHECK_H
#define FFARAWMODULES_LL1CHECK_H

#include "DumpPacket.h"

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class LL1Check : public SubsysReco, public DumpPacket
{
 public:
  LL1Check(const std::string &name = "LL1Check");

  ~LL1Check() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  //  int ResetEvent(PHCompositeNode *topNode) override;

  //  void MyEvtNode(const std::string &name) {m_EvtNodeName = name;}

 private:
  /* std::string m_EvtNodeName = "LL1RAWHIT"; */
  /* std::set<uint64_t> bclk_seen; */
};

#endif  // FFARAWMODULES_LL1CHECK_H
