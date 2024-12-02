#ifndef FUN4ALLRAW_SINGLETRIGGERINPUT_H
#define FUN4ALLRAW_SINGLETRIGGERINPUT_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/InputFileHandler.h>

#include <cstdint>  // for uint64_t
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

class Eventiterator;
class Fun4AllPrdfInputTriggerManager;
class OfflinePacket;
class Packet;
class PHCompositeNode;

class SingleTriggerInput : public Fun4AllBase, public InputFileHandler
{
 public:
  explicit SingleTriggerInput(const std::string &name);
  ~SingleTriggerInput() override;
  virtual Eventiterator *GetEventIterator() { return m_EventIterator; }
  virtual void FillPool(const unsigned int = 1) { return; }
  virtual void RunNumber(const int runno) { m_RunNumber = runno; }
  virtual int RunNumber() const { return m_RunNumber; }
  virtual int fileopen(const std::string &filename) override;
  virtual int fileclose() override;
  virtual int AllDone() const { return m_AllDone; }
  virtual void AllDone(const int i) { m_AllDone = i; }
  virtual int EventNumberOffset(const int packetid);
  virtual void AdjustEventNumberOffset(const int packetid, const int offset);
  virtual void Print(const std::string &what = "ALL") const override;
  virtual void CleanupUsedPackets(const int) { return; }
  virtual bool CheckPoolDepth(const uint64_t bclk);
  virtual void ClearCurrentEvent();
  virtual Eventiterator *GetEventiterator() const { return m_EventIterator; }
  virtual Fun4AllPrdfInputTriggerManager *TriggerInputManager() { return m_TriggerInputMgr; }
  virtual void TriggerInputManager(Fun4AllPrdfInputTriggerManager *in) { m_TriggerInputMgr = in; }
  virtual void CreateDSTNode(PHCompositeNode *) { return; }
  virtual void SubsystemEnum(const int id) { m_SubsystemEnum = id; }
  virtual int SubsystemEnum() const { return m_SubsystemEnum; }
  virtual void ddumppacket(Packet *pkt);
  virtual void enable_ddump(int i = 1) { m_ddump_flag = i; }
  virtual bool ddump_enabled() const { return m_ddump_flag; }
  virtual void DefaultEventNumberOffset(const int i) { m_DefaultEventNumberOffset = i; }
  virtual int AdjustPacketMap(int pktid, int evtoffset);
  virtual bool GetSomeMoreEvents(const unsigned int keep);
  virtual int AdjustEventOffset(int evtoffset);
  virtual void LocalPoolDepth(unsigned int i) { m_LocalPoolDepth = i; }
  virtual unsigned int LocalPoolDepth() const { return m_LocalPoolDepth; }
  virtual void SkipToEvent(const int i) { m_SkipToEvent = i; }
  virtual int SkipToEvent() const { return m_SkipToEvent; }
  virtual void LastEvent(const int i) { m_LastEvent = i; }
  virtual int LastEvent() const { return m_LastEvent; }
  virtual int SetFEMEventRefPacketId(const int pktid);
  virtual int FEMEventRefPacketId() const {return  m_FEMEventRefPacketId;}
  // these ones are used directly by the derived classes, maybe later
  // move to cleaner accessors
 protected:
  std::map<int, std::vector<OfflinePacket *>> m_PacketMap;
  unsigned int m_NumSpecialEvents{0};
  std::set<int> m_EventNumber;
  std::set<int> m_EventStack;

  // we have accessors for these here
 private:
  Eventiterator *m_EventIterator{nullptr};
  Fun4AllPrdfInputTriggerManager *m_TriggerInputMgr{nullptr};
  int m_ddump_flag{0};
  int m_RunNumber{0};
  int m_EventsThisFile{0};
  int m_AllDone{0};
  int m_SubsystemEnum{0};
  int m_DefaultEventNumberOffset{0};
  int m_FEMEventRefPacketId {0};
  int m_SkipToEvent{0};  // we may have negative event numbers but lets not go there right now
  int m_LastEvent{std::numeric_limits<int>::max()};
  unsigned int m_LocalPoolDepth{0};
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
  std::map<int, std::ofstream *> m_PacketDumpFile;
  std::map<int, int> m_PacketDumpCounter;
  std::map<int, int> m_EventNumberOffset;  // packet wise event number offset
};

#endif
