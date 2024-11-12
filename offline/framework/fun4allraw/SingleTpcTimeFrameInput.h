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
class TpcTimeFrameBuilder;

//! Provide TpcTimeFrameBuilder as a unified interface for Fun4AllStreamingInputManager
// NOLINTNEXTLINE(hicpp-special-member-functions)
class SingleTpcTimeFrameInput : public SingleStreamingInput
{
 public:
  explicit SingleTpcTimeFrameInput(const std::string &name);
  ~SingleTpcTimeFrameInput() override;
  void FillPool(const uint64_t minBCO) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  // bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  // bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  void SetBcoRange(const unsigned int i) { m_BcoRange = i; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }

 private:
  const int NTPCPACKETS = 3;

  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};

  //! packet ID -> TimeFrame builder
  std::map<int, TpcTimeFrameBuilder *> m_TpcTimeFrameBuilderMap;
};

#endif
