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

  //! remove used packets matching a given BCO from internal container
  virtual void CleanupUsedPackets(const uint64_t /*BCO*/) {}

  //! remove used packets matching a given BCO from internal container
  /**
   * second parameter is to specify whether BCO has been
   * - succesfully processed or
   * - is dropped
   */
  virtual void CleanupUsedPackets(const uint64_t /*BCO*/ ,bool /*dropped*/) {}

  virtual bool CheckPoolDepth(const uint64_t bclk);
  virtual void ClearCurrentEvent();
  virtual Eventiterator *GetEventiterator() const { return m_EventIterator; }
  virtual Fun4AllStreamingInputManager *StreamingInputManager() { return m_StreamingInputMgr; }
  virtual void StreamingInputManager(Fun4AllStreamingInputManager *in) { m_StreamingInputMgr = in; }
  virtual void CreateDSTNode(PHCompositeNode *) { return; }
  virtual void ConfigureStreamingInputManager() { return; }
  virtual void SubsystemEnum(const int id) { m_SubsystemEnum = id; }
  virtual int SubsystemEnum() const { return m_SubsystemEnum; }
  void MaxBclkDiff(uint64_t ui) { m_MaxBclkSpread = ui; }
  uint64_t MaxBclkDiff() const { return m_MaxBclkSpread; }
  virtual const std::map<int, std::set<uint64_t>> &BclkStackMap() const { return m_BclkStackPacketMap; }
  virtual const std::set<uint64_t> &BclkStack() const { return m_BclkStack; }
  virtual const std::map<uint64_t, std::set<int>> &BeamClockFEE() const { return m_BeamClockFEE; }
  void setHitContainerName(const std::string &name) { m_rawHitContainerName = name; }
  std::string getHitContainerName() const { return m_rawHitContainerName; }
  const std::map<int, std::set<uint64_t>> &getFeeGTML1BCOMap() const { return m_FeeGTML1BCOMap; }

  //! event assembly QA histograms
  virtual void createQAHistos() {}

  //! event assembly QA for a given BCO
  /** TODO: check whether necessary */
  virtual void FillBcoQA(uint64_t /*gtm_bco*/) {};


  void clearPacketBClkStackMap(const int &packetid, const uint64_t& bclk)
  {
    std::set<uint64_t> to_erase;
    auto set = m_BclkStackPacketMap.find(packetid)->second;
      for(auto& bclk_to_erase : set)
      {
        if(bclk_to_erase <= bclk)
        {
          to_erase.insert(bclk_to_erase);
        }
      }
      for(auto& bclk_to_erase : to_erase)
      {
        set.erase(bclk_to_erase);
      }
    }

  void clearFeeGTML1BCOMap(const uint64_t &bclk)
  {
    std::set<uint64_t> toerase;
    for (auto &[key, set] : m_FeeGTML1BCOMap)
    {
      for (auto &ll1bclk : set)
      {
        if (ll1bclk <= bclk)
        {
          // to avoid invalid reads
          toerase.insert(ll1bclk);
        }
      }
      for (auto &bclk_to_erase : toerase)
      {
        set.erase(bclk_to_erase);
      }
    }
  }

 protected:
  std::map<int, std::set<uint64_t>> m_BclkStackPacketMap;
  std::map<int, std::set<uint64_t>> m_FeeGTML1BCOMap;
  std::string m_rawHitContainerName = "";

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
