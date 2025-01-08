#ifndef FUN4ALLRAW_SINGLETPCPOOLINPUT_H
#define FUN4ALLRAW_SINGLETPCPOOLINPUT_H

#include "SingleStreamingInput.h"

#include <array>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class TpcRawHit;
class Packet;

class SingleTpcPoolInput : public SingleStreamingInput
{
 public:
  explicit SingleTpcPoolInput(const std::string &name);
  ~SingleTpcPoolInput() override = default;
  void FillPool(const uint64_t) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents(const uint64_t bclk);

  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  void SetBcoRange(const unsigned int i) { m_BcoRange = i; }
  void SetMaxTpcTimeSamples(const unsigned int i) { m_max_tpc_time_samples = i; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }
  const std::map<int, std::set<uint64_t>> &BclkStackMap() const override { return m_BclkStackPacketMap; }

 private:
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};
  unsigned int m_max_tpc_time_samples{425};
  bool m_skipEarlyEvents{true};
  //! map bco to packet
  std::map<unsigned int, uint64_t> m_packet_bco;

  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<TpcRawHit *>> m_TpcRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif
