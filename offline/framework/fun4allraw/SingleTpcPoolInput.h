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
  ~SingleTpcPoolInput() override;
  void FillPool(const unsigned int) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  void SetBcoRange(const unsigned int i) { m_BcoRange = i; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }

 private:
  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};

  //! map bco to packet
  std::map<unsigned int, uint64_t> m_packet_bco;

  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<TpcRawHit *>> m_TpcRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif
