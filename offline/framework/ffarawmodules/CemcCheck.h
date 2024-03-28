// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_CEMCCHECK_H
#define FFARAWMODULES_CEMCCHECK_H

#include "DumpPacket.h"

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;

class CemcCheck : public SubsysReco, public DumpPacket
{
 public:
  CemcCheck(const std::string &name = "CemcCheck");

  ~CemcCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;


 private:
  /* std::string m_EvtNodeName = "CEMCRAWHIT"; */
  /* std::set<uint64_t> bclk_seen; */
};

#endif  // FFARAWMODULES_CEMCCHECK_H
