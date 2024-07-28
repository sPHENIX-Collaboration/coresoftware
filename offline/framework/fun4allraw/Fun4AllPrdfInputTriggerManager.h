// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLPRDFINPUTTRIGGERMANAGER_H
#define FUN4ALLRAW_FUN4ALLPRDFINPUTTRIGGERMANAGER_H

#include "InputManagerType.h"

#include <fun4all/Fun4AllInputManager.h>

#include <Event/phenixTypes.h>

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

class Event;
class SinglePrdfInput;
class CaloPacket;
class Gl1Packet;
class LL1Packet;
class PHCompositeNode;
class SingleTriggerInput;
class SyncObject;
class OfflinePacket;

class Fun4AllPrdfInputTriggerManager : public Fun4AllInputManager
{
 public:
  Fun4AllPrdfInputTriggerManager(const std::string &name = "DUMMY", const std::string &prdfnodename = "PRDF", const std::string &topnodename = "TOP");
  ~Fun4AllPrdfInputTriggerManager() override;

  int fileopen(const std::string & /* filenam */) override { return 0; }
  // cppcheck-suppress virtualCallInConstructor
  int fileclose() override;
  int run(const int nevents = 0) override;

  void Print(const std::string &what = "ALL") const override;
  int ResetEvent() override;
  int PushBackEvents(const int i) override;
  int GetSyncObject(SyncObject **mastersync) override;
  int SyncIt(const SyncObject *mastersync) override;
  int HasSyncObject() const override { return 1; }
  std::string GetString(const std::string &what) const override;
  void registerTriggerInput(SingleTriggerInput *prdfin, InputManagerType::enu_subsystem system);
  void UpdateEventFoundCounter(const int evtno);
  void UpdateDroppedPacket(const int packetid);
  void AddBeamClock(const int evtno, const int bclk, SinglePrdfInput *prdfin);
  void SetReferenceClock(const int evtno, const int bclk);
  void DitchEvent(const int eventno);
  void ClearAllEvents(const int eventno);
  void SetPoolDepth(unsigned int d) { m_DefaultPoolDepth = d; }
  int FillCemc(const unsigned int nEvents = 2);
  int MoveCemcToNodeTree();
  void AddCemcPacket(int eventno, CaloPacket *pkt);
  int FillGl1(const unsigned int nEvents = 2);
  int MoveGl1ToNodeTree();
  void AddGl1Packet(int eventno, Gl1Packet *gl1pkt);
  int FillLL1(const unsigned int nEvents = 2);
  int MoveLL1ToNodeTree();
  void AddLL1Packet(int eventno, LL1Packet *pkt);
  int FillMbd(const unsigned int nEvents = 2);
  int MoveMbdToNodeTree();
  void AddMbdPacket(int eventno, CaloPacket *mbdpkt);
  int FillHcal(const unsigned int nEvents = 2);
  int MoveHcalToNodeTree();
  void AddHcalPacket(int eventno, CaloPacket *pkt);
  int FillZdc(const unsigned int nEvents = 2);
  int MoveZdcToNodeTree();
  void AddZdcPacket(int eventno, CaloPacket *pkt);
  // the sepd is read together with the zdc in the FillZdc method
  int MoveSEpdToNodeTree();
  void AddSEpdPacket(int eventno, CaloPacket *pkt);
  void InitialPoolDepth(unsigned int n)
  {
    m_InitialPoolDepth = n;
    m_PoolDepth = n;
  }
  void DetermineReferenceEventNumber();
  void ClockDiffFill();
  int ClockDiffCheck();
  void Resync(bool b = true) { m_resync_flag = b; }
  void AddGl1DroppedEvent(int iev) { m_Gl1DroppedEvent.insert(iev); }
  void AddFEMProblemPacket(int i) { m_FEMClockPackets.insert(i); }

 private:
  struct Gl1PacketInfo
  {
    std::map<int, Gl1Packet *> Gl1SinglePacketMap;
    std::map<int, uint64_t> BcoDiffMap;
    unsigned int EventFoundCounter{0};
  };

  struct CaloPacketInfo
  {
    std::map<int, CaloPacket *> CaloSinglePacketMap;
    std::map<int, uint64_t> BcoDiffMap;
    unsigned int EventFoundCounter{0};
  };
  int FillNeedle(std::map<int, CaloPacketInfo>::iterator begin, std::map<int, CaloPacketInfo>::iterator end, const std::string &name = "NONE");
  int ShiftEvents(std::map<int, CaloPacketInfo> &PacketInfoMap, std::map<int, int> &eventoffset, const std::string &name = "NONE");
  int AdjustBcoDiff(std::map<int, CaloPacketInfo> &PacketInfoMap, int packetid, uint64_t bcodiff);
  int DropFirstEvent(std::map<int, CaloPacketInfo> &PacketInfoMap);

  struct LL1PacketInfo
  {
    std::map<int, LL1Packet *> LL1SinglePacketMap;
    std::map<int, uint64_t> BcoDiffMap;
    unsigned int EventFoundCounter{0};
  };
  int FillNeedleLL1(std::map<int, LL1PacketInfo>::iterator begin, std::map<int, LL1PacketInfo>::iterator end, const std::string &name = "NONE");
  int ShiftEventsLL1(std::map<int, LL1PacketInfo> &PacketInfoMap, std::map<int, int> &eventoffset, const std::string &name = "NONE");
  int AdjustBcoDiffLL1(std::map<int, LL1PacketInfo> &PacketInfoMap, int packetid, uint64_t bcodiff);
  int DropFirstEventLL1(std::map<int, LL1PacketInfo> &PacketInfoMap);

  int m_RunNumber{0};
  int m_RefEventNo{std::numeric_limits<int>::min()};
  bool m_gl1_registered_flag{false};
  bool m_mbd_registered_flag{false};
  bool m_cemc_registered_flag{false};
  bool m_hcal_registered_flag{false};
  bool m_ll1_registered_flag{false};
  bool m_zdc_registered_flag{false};
  bool m_resync_flag{false};
  unsigned int m_InitialPoolDepth = 10;
  unsigned int m_DefaultPoolDepth = 10;
  unsigned int m_PoolDepth{m_InitialPoolDepth};
  std::set<int> m_Gl1DroppedEvent;
  std::vector<SingleTriggerInput *> m_TriggerInputVector;
  std::vector<SingleTriggerInput *> m_NoGl1InputVector;
  std::vector<SingleTriggerInput *> m_Gl1InputVector;
  std::vector<SingleTriggerInput *> m_CemcInputVector;
  std::vector<SingleTriggerInput *> m_HcalInputVector;
  std::vector<SingleTriggerInput *> m_LL1InputVector;
  std::vector<SingleTriggerInput *> m_MbdInputVector;
  std::vector<SingleTriggerInput *> m_SEpdInputVector;
  std::vector<SingleTriggerInput *> m_ZdcInputVector;
  SyncObject *m_SyncObject{nullptr};
  PHCompositeNode *m_topNode{nullptr};
  std::map<int, Gl1PacketInfo> m_Gl1PacketMap;
  std::map<int, CaloPacketInfo> m_MbdPacketMap;
  std::map<int, CaloPacketInfo> m_CemcPacketMap;
  std::map<int, CaloPacketInfo> m_HcalPacketMap;
  std::map<int, LL1PacketInfo> m_LL1PacketMap;
  std::map<int, CaloPacketInfo> m_SEpdPacketMap;
  std::map<int, CaloPacketInfo> m_ZdcPacketMap;
  std::map<int, int> m_DroppedPacketMap;
  std::map<int, std::vector<std::pair<int, SinglePrdfInput *>>> m_ClockCounters;
  std::map<int, int> m_RefClockCounters;
  std::vector<uint64_t> m_HayStack;
  std::map<int, std::vector<uint64_t>> m_NeedleMap;
  std::set<int> m_FEMClockPackets;
};

#endif /* FUN4ALL_FUN4ALLPRDFINPUTPOOLMANAGER_H */
