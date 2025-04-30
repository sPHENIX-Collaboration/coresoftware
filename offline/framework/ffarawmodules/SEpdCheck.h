// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_SEPDCHECK_H
#define FFARAWMODULES_SEPDCHECK_H

#include "DumpPacket.h"

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class SEpdCheck : public SubsysReco, public DumpPacket
{
 public:
  SEpdCheck(const std::string &name = "SEpdCheck");

  ~SEpdCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

 private:
  /* std::string m_EvtNodeName = "SEPDRAWHIT"; */
  /* std::set<uint64_t> bclk_seen; */
};

#endif  // FFARAWMODULES_SEPDCHECK_H
