// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_HCALCHECK_H
#define FFARAWMODULES_HCALCHECK_H

#include "DumpPacket.h"

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class HcalCheck : public SubsysReco, public DumpPacket
{
 public:
  HcalCheck(const std::string &name = "HcalCheck");

  ~HcalCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;


 private:
  /* std::string m_EvtNodeName = "HCALRAWHIT"; */
  /* std::set<uint64_t> bclk_seen; */
};

#endif  // FFARAWMODULES_HCALCHECK_H
