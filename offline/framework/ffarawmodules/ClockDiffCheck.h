// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_CLOCKDIFFCHECK_H
#define FFARAWMODULES_CLOCKDIFFCHECK_H

#include "DumpPacket.h"

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class Fun4AllInputManager;
class PHCompositeNode;
class TH1;

class ClockDiffCheck : public SubsysReco, public DumpPacket
{
 public:
  ClockDiffCheck(const std::string &name = "ClockDiffCheck");

  ~ClockDiffCheck() override {}

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;


 private:
  std::map<unsigned int, std::tuple<unsigned int, TH1 *>> m_PacketStuffMap;
  /* std::string m_EvtNodeName = "CLOCKDIFFRAWHIT"; */
  /* std::set<uint64_t> bclk_seen; */
};

#endif  // FFARAWMODULES_CLOCKDIFFCHECK_H
