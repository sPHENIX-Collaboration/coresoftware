#ifndef FUN4ALLRAW_SINGLETRIGGERINPUT_H
#define FUN4ALLRAW_SINGLETRIGGERINPUT_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/InputFileHandler.h>

#include <cstdint>  // for uint64_t
#include <map>
#include <set>
#include <string>

class Eventiterator;
class Fun4AllStreamingInputManager;
class Fun4AllPrdfInputTriggerManager;
class PHCompositeNode;

class SingleTriggerInput : public Fun4AllBase, public InputFileHandler
{
 public:
  explicit SingleTriggerInput(const std::string &name);
  ~SingleTriggerInput() override;
  virtual Eventiterator *GetEventIterator() { return m_EventIterator; }
  virtual void FillPool(const uint64_t) { return; }
  virtual void FillPool(const unsigned int = 1) { return; }
  virtual void RunNumber(const int runno) { m_RunNumber = runno; }
  virtual int RunNumber() const { return m_RunNumber; }
  virtual int fileopen(const std::string &filename) override;
  virtual int fileclose() override;
  virtual int AllDone() const { return m_AllDone; }
  virtual void AllDone(const int i) { m_AllDone = i; }
  virtual void EventNumberOffset(const int i) { m_EventNumberOffset = i; }
  virtual void Print(const std::string &what = "ALL") const override;
  virtual void CleanupUsedPackets(const int) { return; }
  virtual bool CheckPoolDepth(const uint64_t bclk);
  virtual void ClearCurrentEvent();
  virtual Eventiterator *GetEventiterator() const { return m_EventIterator; }
  virtual Fun4AllStreamingInputManager *StreamingInputManager() { return m_StreamingInputMgr; }
  virtual void StreamingInputManager(Fun4AllStreamingInputManager *in) { m_StreamingInputMgr = in; }
  virtual Fun4AllPrdfInputTriggerManager *TriggerInputManager() { return m_TriggerInputMgr; }
  virtual void TriggerInputManager(Fun4AllPrdfInputTriggerManager *in) { m_TriggerInputMgr = in; }
  virtual void CreateDSTNode(PHCompositeNode *) { return; }
  virtual void ConfigureStreamingInputManager() { return; }
  virtual void SubsystemEnum(const int id) { m_SubsystemEnum = id; }
  virtual int SubsystemEnum() const { return m_SubsystemEnum; }

 private:
  Eventiterator *m_EventIterator = nullptr;
  Fun4AllStreamingInputManager *m_StreamingInputMgr = nullptr;
  Fun4AllPrdfInputTriggerManager *m_TriggerInputMgr{nullptr};
  unsigned int m_EventNumberOffset = 1;  // packet event counters start at 0 but we start with event number 1
  int m_RunNumber = 0;
  int m_EventsThisFile = 0;
  int m_AllDone = 0;
  int m_SubsystemEnum{0};
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif
