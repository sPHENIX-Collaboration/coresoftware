#ifndef FUN4ALLRAW_SINGLEINTTPOOLINPUT_H
#define FUN4ALLRAW_SINGLEINTTPOOLINPUT_H

#include "SingleStreamingInput.h"

#include <array>
#include <cstdint>  // for uint64_t
#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>
class InttRawHit;
class Packet;
class PHCompositeNode;
class intt_pool;

class SingleInttPoolInput : public SingleStreamingInput
{
 public:
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
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }
  const std::set<uint64_t> &BclkStack() const override { return m_BclkStack; }
  const std::map<uint64_t, std::set<int>> &BeamClockFEE() const override { return m_BeamClockFEE; }
  void streamingMode(const bool isStreaming)
  {
    if(isStreaming)
    {
      SetNegativeBco(120 - 23);
      SetBcoRange(500);
      std::cout << "INTT set to streaming event combining"<<std::endl;
      return;
    }

    SetNegativeBco(1);
    SetBcoRange(2);
    std::cout << "INTT set to triggered event combining" << std::endl;
  }

 private:
  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};
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
