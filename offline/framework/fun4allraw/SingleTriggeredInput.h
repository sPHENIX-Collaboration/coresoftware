#ifndef FUN4ALLRAW_SINGLETRIGGEREDINPUT_H
#define FUN4ALLRAW_SINGLETRIGGEREDINPUT_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/InputFileHandler.h>

#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>
#include <Event/packet.h>

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
  /**
 * Set the current run number used for event processing.
 * @param runno The run number to assign; this value will be associated with subsequently read events.
 */
virtual void RunNumber(const int runno) { m_RunNumber = runno; }
  virtual int RunNumber() const { return m_RunNumber; }
  virtual int fileopen(const std::string &filename) override;
  /**
   * Close the currently opened input file.
   * @returns 0 on success, non-zero on error.
   */
  
  /**
   * Report whether processing is complete.
   * @returns The completion flag: nonzero if processing is complete, zero otherwise.
   */
  
  /**
   * Set the completion flag for processing.
   * @param i Completion flag value; nonzero indicates processing is complete.
   */
  virtual int fileclose() override;
  virtual int AllDone() const { return m_AllDone; }
  virtual void AllDone(const int i) { m_AllDone = i; }
  virtual int FilesDone() const { return m_FilesDone; }
  /**
 * Set the count of files that have been completed.
 * @param i Number of completed files to store as the current files-done count.
 */
virtual void FilesDone(const int i) { m_FilesDone = i; }
  /**
 * Set the event alignment problem flag.
 *
 * @param i Flag value where non-zero indicates an event alignment problem, zero clears it.
 */
virtual void EventAlignmentProblem(const int i) { m_EventAlignmentProblem = i; }
  virtual int EventAlignmentProblem() const { return m_EventAlignmentProblem; }
  /**
 * Set the current event number used by this input instance.
 * @param i Event number to set; updates the internal event counter state.
 */
virtual void EventNumber(const int i) { m_EventNumber = i; }
  virtual int EventNumber() const { return m_EventNumber; }
  virtual void CreateDSTNodes(Event *evt);
  // these ones are used directly by the derived classes, maybe later
  // move to cleaner accessors
  virtual int FillEventVector();
  virtual void FillPacketClock(Event *evt, Packet *pkt, size_t event_index);
  virtual int ReadEvent();
  virtual SingleTriggeredInput *Gl1Input() { return m_Gl1Input; }
  virtual void Gl1Input(SingleTriggeredInput *input) { m_Gl1Input = input; }
  virtual uint64_t GetClock(Event *evt, int pid);
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
  virtual bool CheckFemDiffIdx(int pid, size_t index, const std::deque<Event*>& events, uint64_t gl1diffidx);
  virtual bool CheckPoolAlignment(int pid, const std::array<uint64_t, pooldepth>& sebdiff, const std::array<uint64_t, pooldepth>& gl1diff, std::vector<int>& bad_indices, int& shift, bool& CurrentPoolLastDiffBad, bool PrevPoolLastDiffBad);
  virtual bool FemClockAlignment(int pid, const std::deque<Event*>& events, const std::array<uint64_t, pooldepth>& gl1diff);

 protected:
  PHCompositeNode *m_topNode{nullptr};
  // lined up like this:
  // Event * | previous event beam clock | clock diff to previous event
  // keeping previous beam clock just eases the looping, we want to be able to have
  // the accompanying diff to the previous beam clock with this event, so any mismatch
  // gives us the event index in the deque which is off
  std::deque<Event *> m_EventDeque;
  std::map<int, std::deque<Event*>> m_PacketEventDeque;
  std::map<int, Event*> m_PacketEventBackup;
  std::map<int, int> m_PacketShiftOffset;
  std::array<uint64_t, pooldepth + 1> m_bclkarray{};  // keep the last bco from previous loop
  std::array<uint64_t, pooldepth> m_bclkdiffarray{};
  std::map<int, std::array<uint64_t, pooldepth + 1>> m_bclkarray_map;
  std::map<int, std::array<uint64_t, pooldepth>>     m_bclkdiffarray_map;
  std::set<int> m_PacketSet;
  static uint64_t ComputeClockDiff(uint64_t curr, uint64_t prev) { return (curr - prev) & 0xFFFFFFFF; }


 private:
  Eventiterator *m_EventIterator{nullptr};
  SingleTriggeredInput *m_Gl1Input{nullptr};
  /**
 * Flag indicating whether processing is complete.
 *
 * When non-zero, the input has finished processing and no further events or files will be read.
 */
int m_AllDone{0};
  /**
 * @brief Current event sequence number within the input stream.
 *
 * Tracks the number of events processed by this SingleTriggeredInput instance;
 * used for bookkeeping, diagnostics, and run/event lifecycle management.
 */
uint64_t m_Event{0};
  /**
 * Current event number within the input stream.
 *
 * Stores the index of the most recently processed or loaded event and is used to track progression through input events.
 */
int m_EventNumber{0};
  /**
 * Indicates whether an event alignment problem has been detected.
 *
 * 0 means no alignment problem; a non-zero value indicates an alignment problem.
 */
int m_EventAlignmentProblem{0};
  /**
 * Number of input files that have been completely processed.
 *
 * Tracks how many files the input handler has finished handling; used to report progress and determine completion.
 */
int m_FilesDone{0};
  int m_LastEvent{std::numeric_limits<int>::max()};
  int m_ProblemEvent{-1};
  int m_RepresPacket{-1};
  int m_RunNumber{0};
  int m_max_alignment_retries{5};
  bool firstcall{true};
  bool firstclockcheck{true};
  bool m_KeepPacketsFlag{false};
  bool m_packetclk_copy_runs{false};
  std::set<int> m_CorrectCopiedClockPackets;
  std::map<int, std::set<int>> m_DitchPackets;
  std::set<int> m_FEMEventNrSet;
  std::set<int> m_OverrideWithRepClock;
  std::map<int, int> m_PacketAlignmentFailCount;
  std::map<int, bool> m_PacketAlignmentProblem;
  std::map<int, bool> m_PrevPoolLastDiffBad;
  std::map<int, uint64_t> m_PreviousValidBCOMap;
  /**
 * Monotonic counter of processed events for this SingleTriggeredInput instance.
 *
 * Tracks the total number of events observed/handled by the object across its
 * lifetime; incremented as events are read and used for bookkeeping and
 * diagnostics.
 */
long long eventcounter{0};
};

#endif