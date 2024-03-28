// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_ZDCCHECK_H
#define FFARAWMODULES_ZDCCHECK_H

#include "DumpPacket.h"

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class ZdcCheck : public SubsysReco, public DumpPacket
{
 public:
  ZdcCheck(const std::string &name = "ZdcCheck");

  ~ZdcCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;


 private:
  /* std::string m_EvtNodeName = "ZDCRAWHIT"; */
  /* std::set<uint64_t> bclk_seen; */
};

#endif  // FFARAWMODULES_ZDCCHECK_H
