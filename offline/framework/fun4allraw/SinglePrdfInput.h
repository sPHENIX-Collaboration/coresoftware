#ifndef FUN4ALLRAW_SINGLEPRDFINPUT_H
#define FUN4ALLRAW_SINGLEPRDFINPUT_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/InputFileHandler.h>

#include <array>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>  // for pair
#include <vector>

class Eventiterator;
class Fun4AllPrdfInputPoolManager;
class Fun4AllPrdfInputTriggerManager;
class Packet;

class SinglePrdfInput : public Fun4AllBase, public InputFileHandler
{
 public:
  explicit SinglePrdfInput(const std::string &name, Fun4AllPrdfInputPoolManager *inman);
  explicit SinglePrdfInput(const std::string &name, Fun4AllPrdfInputTriggerManager *inman);
  ~SinglePrdfInput() override;
  Eventiterator *GetEventIterator() { return m_EventIterator; }
  virtual void FillPool(const unsigned int nevents);
  int RunNumber() const { return m_RunNumber; }
  void RunNumber(const int irun) { m_RunNumber = irun; }
  int fileopen(const std::string &filename) override;
  int fileclose() override;
  int AllDone() const { return m_AllDone; }
  void AllDone(const int i) { m_AllDone = i; }
  void EventNumberOffset(const int i) { m_EventNumberOffset = i; }  // if beam clk are out of sync, tweak this one
  int EventNumberOffset() const { return m_EventNumberOffset; }
  Fun4AllPrdfInputPoolManager *InputMgr() { return m_InputMgr; }
  Fun4AllPrdfInputTriggerManager *TriggerInputMgr() { return m_TriggerInputMgr; }
  void MakeReference(const bool b);
  bool ReferenceFlag() const { return m_MeReferenceFlag; }

 private:
  int majority_eventnumber();
  int majority_beamclock();
  void adjust_eventnumber_offset(const int decided_evtno);
  struct PacketInfo
  {
    std::vector<Packet *> PacketVector;
    unsigned int EventFoundCounter = 0;
  };
  Eventiterator *m_EventIterator = nullptr;
  Fun4AllPrdfInputPoolManager *m_InputMgr = nullptr;
  Fun4AllPrdfInputTriggerManager *m_TriggerInputMgr = nullptr;
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;
  int m_EventNumberOffset = 1;               // packet event counters start at 0 but we start with event number 1
  int *m_PacketEventNumberOffset = nullptr;  // packet event counters start at 0 but we start with event number 1
  int m_RunNumber = 0;
  int m_EventsThisFile = 0;
  int m_AllDone = 0;
  bool m_MeReferenceFlag = false;
  std::map<int, std::vector<Packet *>> m_PacketMap;
  std::set<int> m_EvtSet;
  std::vector<std::pair<int, int>> m_Event;
  std::array<unsigned int, 100> rollover{};
  std::array<int, 100> previous_eventnumber{};
};

#endif
