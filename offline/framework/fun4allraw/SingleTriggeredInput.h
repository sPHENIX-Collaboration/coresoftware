#ifndef FUN4ALLRAW_SINGLETRIGGEREDINPUT_H
#define FUN4ALLRAW_SINGLETRIGGEREDINPUT_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/InputFileHandler.h>

#include <array>
#include <cstdint>  // for uint64_t
#include <deque>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

class Event;
class Eventiterator;
class OfflinePacket;
class PHCompositeNode;

class SingleTriggeredInput : public Fun4AllBase, public InputFileHandler
{
 public:
  static constexpr size_t pooldepth{10};  // number of events which are read in in one go
  explicit SingleTriggeredInput(const std::string &name);
  ~SingleTriggeredInput() override;
  virtual Eventiterator *GetEventIterator() { return m_EventIterator; }
  virtual void FillPool();
  virtual void RunNumber(const int runno) { m_RunNumber = runno; }
  virtual int RunNumber() const { return m_RunNumber; }
  virtual void EventNumber(const int i) { m_EventNumber = i; }
  virtual int EventNumber() const { return m_EventNumber; }
  virtual int EventsInThisFile() const { return m_EventsThisFile; }
  virtual int fileopen(const std::string &filename) override;
  virtual int fileclose() override;
  virtual int AllDone() const { return m_AllDone; }
  virtual void AllDone(const int i) { m_AllDone = i; }
  virtual int FilesDone() const { return m_FilesDone; }
  virtual void FilesDone(const int i) { m_FilesDone = i; }
  virtual void EventAlignmentProblem(const int i) { m_EventAlignmentProblem = i; }
  virtual int EventAlignmentProblem() const { return m_EventAlignmentProblem; }
  virtual void CreateDSTNodes(Event *evt);
  // these ones are used directly by the derived classes, maybe later
  // move to cleaner accessors
  virtual int FillEventVector();
  virtual int ReadEvent();
  virtual SingleTriggeredInput *Gl1Input() { return m_Gl1Input; }
  virtual void Gl1Input(SingleTriggeredInput *input) { m_Gl1Input = input; }
  virtual uint64_t GetClock(Event *evt);
  virtual std::array<uint64_t, pooldepth>::const_iterator clkdiffbegin() { return m_bclkdiffarray.begin(); }
  virtual std::array<uint64_t, pooldepth>::const_iterator clkdiffend() { return m_bclkdiffarray.end(); }
  virtual std::array<uint64_t, pooldepth>::const_iterator beginclock() { return m_bclkarray.begin(); }
  virtual void KeepPackets() { m_KeepPacketsFlag = true; }
  virtual bool KeepMyPackets() const { return m_KeepPacketsFlag; }
  void topNode(PHCompositeNode *topNode) { m_topNode = topNode; }
  PHCompositeNode *topNode() { return m_topNode; }
  virtual void FakeProblemEvent(const int ievent) { m_ProblemEvent = ievent; }
  virtual int FemEventNrClockCheck(OfflinePacket *calopkt);
  void dumpdeque();
  int checkfirstsebevent();

 protected:
  PHCompositeNode *m_topNode{nullptr};
  // lined up like this:
  // Event * | previous event beam clock | clock diff to previous event
  // keeping previous beam clock just eases the looping, we want to be able to have
  // the accompanying diff to the previous beam clock with this event, so any mismatch
  // gives us the event index in the deque which is off
  std::deque<Event *> m_EventDeque;
  std::array<uint64_t, pooldepth + 1> m_bclkarray{};  // keep the last bco from previous loop
  std::array<uint64_t, pooldepth> m_bclkdiffarray{};
  std::set<int> m_PacketSet;

 private:
  Eventiterator *m_EventIterator{nullptr};
  SingleTriggeredInput *m_Gl1Input{nullptr};
  uint64_t m_Event{0};
  int m_RunNumber{0};
  int m_EventNumber{0};
  int m_EventsThisFile{0};
  int m_AllDone{0};
  int m_FilesDone{0};
  int m_EventAlignmentProblem{0};
  int m_ProblemEvent{-1};
  int m_LastEvent{std::numeric_limits<int>::max()};
  bool firstcall{true};
  bool m_KeepPacketsFlag{false};
  bool m_DitchPackets{false};
  std::set<int> m_FEMEventNrSet;
};

#endif
