#ifndef FUN4ALLRAW_SINGLEEVTINPUT_H
#define FUN4ALLRAW_SINGLEEVTINPUT_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/InputFileHandler.h>

#include <array>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class Eventiterator;
class Fun4AllEvtInputPoolManager;
class Packet;

class SingleEvtInput : public Fun4AllBase, public InputFileHandler
{
 public:
  explicit SingleEvtInput(const std::string &name, Fun4AllEvtInputPoolManager *inman);
  ~SingleEvtInput() override;
  Eventiterator *GetEventIterator() { return m_EventIterator; }
  void FillPool(const unsigned int nevents = 1);
  int RunNumber() const { return m_RunNumber; }
  int fileopen(const std::string &filename) override;
  int fileclose() override;
  int AllDone() const { return m_AllDone; }
  void AllDone(const int i) { m_AllDone = i; }
  void EventNumberOffset(const int i) { m_EventNumberOffset = i; }
  void Print(const std::string &what = "ALL") const override;
  void CleanupUsedPackets(const uint64_t bclk);
  bool CheckPoolDepth(const uint64_t bclk);
  void ClearCurrentEvent();

 private:
  Eventiterator *m_EventIterator = nullptr;
  Fun4AllEvtInputPoolManager *m_InputMgr = nullptr;
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;
  unsigned int m_EventNumberOffset = 1;  // packet event counters start at 0 but we start with event number 1
  int m_RunNumber = 0;
  int m_EventsThisFile = 0;
  int m_AllDone = 0;
  std::array<uint64_t, 14> m_PreviousClock{};
  std::array<uint64_t, 14> m_Rollover{};
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<Packet *>> m_PacketStorageMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif
