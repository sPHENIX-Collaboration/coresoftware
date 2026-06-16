#ifndef FUN4ALLRAW_SingleTpcTimeFrameInput_H
#define FUN4ALLRAW_SingleTpcTimeFrameInput_H

#include "SingleStreamingInput.h"

#include <array>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class TpcRawHit;
class Packet;
class TpcTimeFrameBuilderBase;
class PHTimer;
class TH1;
class TH2;

//! Provide TpcTimeFrameBuilder as a unified interface for Fun4AllStreamingInputManager
// NOLINTNEXTLINE(hicpp-special-member-functions)
class SingleTpcTimeFrameInput : public SingleStreamingInput
{
 public:
  explicit SingleTpcTimeFrameInput(const std::string &name);
  ~SingleTpcTimeFrameInput() override;
  void FillPool(const uint64_t targetBCO) override;
  int FillPoolStatus() const override { return m_FillPoolStatus; }
  void CleanupUsedPackets(const uint64_t bclk) override;
  // bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  // bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  void SetBcoRange(const unsigned int i) { m_BcoRange = i; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }

  void AddPacketID(const int packetID) { m_SelectedPacketIDs.insert(packetID); }

  void setDigitalCurrentDebugTTreeName(const std::string &name)
  {
    m_digitalCurrentDebugTTreeName = name;
  }

  void setBXCounterSyncCDBTTreeName(const std::string &name)
  {
    m_bxCounterSyncCDBTTreeName = name;
  }

 private:
  const int NTPCPACKETS = 3;

  // in BCO, limit caching to a quarter of FEE clock rollover or 7ms, to avoid memory over usage when trigger jumped by a long time
  static constexpr uint64_t kUsedPacketsCachingLimit = (1<<20)/4/4;  

  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};

  //! packet ID -> TimeFrame builder
  std::map<int, TpcTimeFrameBuilderBase *> m_TpcTimeFrameBuilderMap;
  std::map<int, int> m_TpcTimeFrameBuilderHitFormatMap;
  std::set<int> m_SelectedPacketIDs;

  TH1 *m_hNorm = nullptr;

  PHTimer *m_FillPoolTimer = nullptr;
  PHTimer *m_getNextEventTimer = nullptr;
  PHTimer *m_ProcessPacketTimer = nullptr;
  PHTimer *m_getTimeFrameTimer = nullptr;

  // NOLINTNEXTLINE(hicpp-special-member-functions)
  class TimeTracker
  {
   public:
    TimeTracker(PHTimer *timer, const std::string &name, TH1 *hout);
    virtual ~TimeTracker();
    void stop();

   private:
    PHTimer *m_timer = nullptr;
    std::string m_name;
    TH1 *m_hNorm = nullptr;
    bool stopped = false;
  };

  int m_FillPoolStatus{0};
  std::string m_digitalCurrentDebugTTreeName;
  std::string m_bxCounterSyncCDBTTreeName;
};

#endif
