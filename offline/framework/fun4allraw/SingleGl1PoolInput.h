#ifndef FUN4ALLRAW_SINGLEGL1POOLINPUT_H
#define FUN4ALLRAW_SINGLEGL1POOLINPUT_H

#include "SingleStreamingInput.h"

#include <cstdint>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class Gl1RawHit;
class PHCompositeNode;

class SingleGl1PoolInput : public SingleStreamingInput
{
 public:
  explicit SingleGl1PoolInput(const std::string &name);
  ~SingleGl1PoolInput() override;
  void FillPool(const unsigned int) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  void SetBcoRange(const unsigned int i) { m_BcoRange = i; }
  //  void ConfigureStreamingInputManager() override;

 private:
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};

  //! map bco to packet
  std::map<unsigned int, uint64_t> m_packet_bco;

  std::set<uint64_t> m_BeamClockFEE;
  std::map<uint64_t, std::vector<Gl1RawHit *>> m_Gl1RawHitMap;
  std::set<uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif
