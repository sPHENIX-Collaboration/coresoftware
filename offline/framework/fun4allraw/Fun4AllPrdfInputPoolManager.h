// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALLRAW_FUN4ALLPRDFINPUTPOOLMANAGER_H
#define FUN4ALLRAW_FUN4ALLPRDFINPUTPOOLMANAGER_H

#include <fun4all/Fun4AllInputManager.h>

#include <Event/phenixTypes.h>

#include <map>
#include <string>

class Event;
class SinglePrdfInput;
class oEvent;
class Packet;
class PHCompositeNode;
class SyncObject;

class Fun4AllPrdfInputPoolManager : public Fun4AllInputManager
{
 public:
  Fun4AllPrdfInputPoolManager(const std::string &name = "DUMMY", const std::string &prdfnodename = "PRDF", const std::string &topnodename = "TOP");
  ~Fun4AllPrdfInputPoolManager() override;
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
  void AddPacket(const int evtno, Packet *p);
  void UpdateEventFoundCounter(const int evtno);
  void UpdateDroppedPacket(const int packetid);
  void AddBeamClock(const int evtno, const int bclk, SinglePrdfInput *prdfin);
  void SetReferenceClock(const int evtno, const int bclk);
  void SetReferenceInputMgr(SinglePrdfInput *inp) { m_RefPrdfInput = inp; }
  void CreateBclkOffsets();
  int CalcDiffBclk(const int bclk1, const int bclk2);
  void DitchEvent(const int eventno);
  void Resynchronize();
  void ClearAllEvents();
  void SetPoolDepth(unsigned int d) {m_PoolDepth = d;}

 private:
  struct PacketInfo
  {
    std::vector<Packet *> PacketVector;
    unsigned int EventFoundCounter = 0;
  };

  struct SinglePrdfInputInfo
  {
    int bclkoffset = 0;
  };

  bool m_StartUpFlag = true;
  int m_RunNumber = 0;
  unsigned int m_PoolDepth = 100;
  unsigned int m_InitialPoolDepth = 20;
  std::vector<SinglePrdfInput *> m_PrdfInputVector;
  SyncObject *m_SyncObject = nullptr;
  PHCompositeNode *m_topNode = nullptr;
  Event *m_Event = nullptr;
  PHDWORD workmem[4 * 1024 * 1024] = {};
  oEvent *oph = nullptr;
  SinglePrdfInput *m_RefPrdfInput = nullptr;
  std::map<int, PacketInfo> m_PacketMap;
  std::string m_PrdfNodeName;
  std::map<int, int> m_DroppedPacketMap;
  std::map<int, std::vector<std::pair<int, SinglePrdfInput *>>> m_ClockCounters;
  std::map<int, int> m_RefClockCounters;
  std::map<SinglePrdfInput *, SinglePrdfInputInfo> m_SinglePrdfInputInfo;
};

#endif /* FUN4ALL_FUN4ALLPRDFINPUTPOOLMANAGER_H */
