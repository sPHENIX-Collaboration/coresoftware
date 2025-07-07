// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FFARAWMODULES_CLOCKDIFFCHECK_H
#define FFARAWMODULES_CLOCKDIFFCHECK_H

#include "DumpPacket.h"

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>
#include <vector>

class CaloPacket;
class CaloPacketContainer;
class Fun4AllInputManager;
class PHCompositeNode;
class TH1;

class ClockDiffCheck : public SubsysReco, public DumpPacket
{
 public:
  ClockDiffCheck(const std::string &name = "ClockDiffCheck");

  ~ClockDiffCheck() override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  void FillCaloClockDiff(CaloPacketContainer *pktcont);
  void FillCaloClockDiffSngl(CaloPacket *calopkt);
  void FillPacketDiff(OfflinePacket *pkt);

  static bool CheckFemEventNr(CaloPacket *calopkt);
  void set_delBadPkts(bool newDelBadPkts)
  {
    delBadPkts = newDelBadPkts;
  }

  bool get_delBadPkts()
  {
    return delBadPkts;
  }
  int count{0};

 private:
  bool delBadPkts{false};
  std::map<unsigned int, std::tuple<uint64_t, uint64_t, uint64_t, TH1 *, bool>> m_PacketStuffMap;
  std::vector<std::string> m_PacketNodeNames;
};

#endif  // FFARAWMODULES_CLOCKDIFFCHECK_H
