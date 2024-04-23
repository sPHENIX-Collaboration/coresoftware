// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_MBDCHECK_H
#define FFARAWMODULES_MBDCHECK_H

#include "DumpPacket.h"

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class MbdCheck : public SubsysReco, public DumpPacket
{
 public:
  MbdCheck(const std::string &name = "MbdCheck");

  ~MbdCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  //  int ResetEvent(PHCompositeNode *topNode) override;

  //  void MyEvtNode(const std::string &name) {m_EvtNodeName = name;}

 private:
  /* std::string m_EvtNodeName = "MBDRAWHIT"; */
  /* std::set<uint64_t> bclk_seen; */
};

#endif  // FFARAWMODULES_MBDCHECK_H
