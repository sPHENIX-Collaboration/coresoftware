// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLPRDFINPUTTRIGGERMANAGER_H
#define FUN4ALLRAW_FUN4ALLPRDFINPUTTRIGGERMANAGER_H

#include "InputManagerType.h"

#include <fun4all/Fun4AllInputManager.h>

#include <Event/phenixTypes.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

class Event;
class SinglePrdfInput;
class oEvent;
class Packet;
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
  SinglePrdfInput *AddPrdfInputList(const std::string &listfile);
  SinglePrdfInput *AddPrdfInputFile(const std::string &filename);
  SinglePrdfInput *registerPrdfInput(SinglePrdfInput *prdfin);
  void registerTriggerInput(SingleTriggerInput *prdfin, InputManagerType::enu_subsystem system);
  void AddPacket(const int evtno, Packet *p);
  void UpdateEventFoundCounter(const int evtno);
  void UpdateDroppedPacket(const int packetid);
  void AddBeamClock(const int evtno, const int bclk, SinglePrdfInput *prdfin);
  void SetReferenceClock(const int evtno, const int bclk);
  void SetReferenceInputMgr(SinglePrdfInput *inp) { m_RefPrdfInput = inp; }
  void CreateBclkOffsets();
  uint64_t CalcDiffBclk(const uint64_t bclk1, const uint64_t bclk2);
  void DitchEvent(const int eventno);
  void Resynchronize();
  void ClearAllEvents();
  void SetPoolDepth(unsigned int d) { m_PoolDepth = d; }
  int FillGl1();
  void AddGl1Packet(int eventno, OfflinePacket *gl1pkt);
  int FillMbd();
  void AddMbdPacket(int eventno, OfflinePacket *mbdpkt);

 private:
  struct PacketInfo
  {
    std::vector<Packet *> PacketVector;
    unsigned int EventFoundCounter = 0;
  };

  struct SinglePrdfInputInfo
  {
    uint64_t bclkoffset = 0;
  };

  struct Gl1PacketInfo
  {
    std::vector<OfflinePacket *> Gl1PacketVector;
    unsigned int EventFoundCounter = 0;
  };

  struct MbdPacketInfo
  {
    std::vector<OfflinePacket *> MbdPacketVector;
    unsigned int EventFoundCounter = 0;
  };

  bool m_StartUpFlag = true;
  int m_RunNumber{0};
  bool m_gl1_registered_flag{false};
  bool m_mbd_registered_flag{false};
  unsigned int m_PoolDepth = 100;
  unsigned int m_InitialPoolDepth = 20;
  int m_RefEventNo{0};
  std::vector<SinglePrdfInput *> m_PrdfInputVector;
  std::vector<SingleTriggerInput *> m_TriggerInputVector;
  std::vector<SingleTriggerInput *> m_Gl1InputVector;
  std::vector<SingleTriggerInput *> m_MbdInputVector;
  SyncObject *m_SyncObject = nullptr;
  PHCompositeNode *m_topNode = nullptr;
  Event *m_Event = nullptr;
  PHDWORD workmem[4 * 1024 * 1024] = {};
  oEvent *oph = nullptr;
  SinglePrdfInput *m_RefPrdfInput = nullptr;
  std::map<int, PacketInfo> m_PacketMap;
  std::map<int, Gl1PacketInfo> m_Gl1PacketMap;
  std::map<int, MbdPacketInfo> m_MbdPacketMap;
  std::string m_PrdfNodeName;
  std::map<int, int> m_DroppedPacketMap;
  std::map<int, std::vector<std::pair<int, SinglePrdfInput *>>> m_ClockCounters;
  std::map<int, int> m_RefClockCounters;
  std::map<SinglePrdfInput *, SinglePrdfInputInfo> m_SinglePrdfInputInfo;
};

#endif /* FUN4ALL_FUN4ALLPRDFINPUTPOOLMANAGER_H */
