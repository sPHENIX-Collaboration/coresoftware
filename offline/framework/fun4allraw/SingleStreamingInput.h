#ifndef FUN4ALLRAW_SINGLESTREAMINGINPUT_H
#define FUN4ALLRAW_SINGLESTREAMINGINPUT_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/InputFileHandler.h>

#include <cstdint>  // for uint64_t
#include <map>
#include <set>
#include <string>

class Eventiterator;
class Fun4AllEvtInputPoolManager;
class Fun4AllStreamingInputManager;
class PHCompositeNode;

class SingleStreamingInput : public Fun4AllBase, public InputFileHandler
{
 public:
  explicit SingleStreamingInput(const std::string &name);
  ~SingleStreamingInput() override;
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
  virtual void CleanupUsedPackets(const uint64_t) { return; }
  virtual bool CheckPoolDepth(const uint64_t bclk);
  virtual void ClearCurrentEvent();
  virtual Eventiterator *GetEventiterator() const { return m_EventIterator; }
  /* virtual Fun4AllEvtInputPoolManager *InputManager() { return m_InputMgr; } */
  /* virtual void InputManager(Fun4AllEvtInputPoolManager *in) { m_InputMgr = in; } */
  virtual Fun4AllStreamingInputManager *StreamingInputManager() { return m_StreamingInputMgr; }
  virtual void StreamingInputManager(Fun4AllStreamingInputManager *in) { m_StreamingInputMgr = in; }
  virtual void CreateDSTNode(PHCompositeNode *) { return; }
  virtual void ConfigureStreamingInputManager() { return; }
  virtual void SubsystemEnum(const int id) { m_SubsystemEnum = id; }
  virtual int SubsystemEnum() const { return m_SubsystemEnum; }
  void MaxBclkDiff(uint64_t ui) { m_MaxBclkSpread = ui; }
  uint64_t MaxBclkDiff() const { return m_MaxBclkSpread; }

 private:
  Eventiterator *m_EventIterator{nullptr};
  //  Fun4AllEvtInputPoolManager *m_InputMgr {nullptr};
  Fun4AllStreamingInputManager *m_StreamingInputMgr{nullptr};
  uint64_t m_MaxBclkSpread{1000000};
  unsigned int m_EventNumberOffset{1};  // packet event counters start at 0 but we start with event number 1
  int m_RunNumber{0};
  int m_EventsThisFile{0};
  int m_AllDone{0};
  int m_SubsystemEnum{0};
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif
