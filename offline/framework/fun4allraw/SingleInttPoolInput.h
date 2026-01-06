#ifndef FUN4ALLRAW_SINGLEINTTPOOLINPUT_H
#define FUN4ALLRAW_SINGLEINTTPOOLINPUT_H

#include "SingleStreamingInput.h"

#include <array>
#include <cstdint>  // for uint64_t
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

class InttRawHit;
class Packet;
class PHCompositeNode;
class intt_pool;

class SingleInttPoolInput : public SingleStreamingInput
{
 public:
  enum InttStreamingMode
  {
    UNDEFINED = 0,
    STREAMING = 1,
    TRIGGERED = -1
  };
  explicit SingleInttPoolInput(const std::string &name);
  ~SingleInttPoolInput() override;
  void FillPool(const uint64_t minBCO) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents(const uint64_t ibclk);
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

  void SetBcoRange(const unsigned int value) { m_BcoRange = value; }
  unsigned int GetBcoRange() const { return m_BcoRange; }
  void ConfigureStreamingInputManager() override { return; }
  void ConfigureStreamingInputManagerLocal(const int runnumber);
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }
  unsigned int GetNegativeBco() const { return m_NegativeBco; }
  const std::set<uint64_t> &BclkStack() const override { return m_BclkStack; }
  const std::map<uint64_t, std::set<int>> &BeamClockFEE() const override { return m_BeamClockFEE; }

  void streamingMode(const bool isStreaming);
  bool IsStreaming(int runnumber);

 private:
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};
  int m_SavedRunNumber{0};
  int m_StreamingFlag{InttStreamingMode::UNDEFINED};
  bool m_SkipEarlyEvents{true};
  std::array<uint64_t, 14> m_PreviousClock{};
  std::array<uint64_t, 14> m_Rollover{};
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<InttRawHit *>> m_InttRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;

  std::map<int, intt_pool *> poolmap;
};

#endif
