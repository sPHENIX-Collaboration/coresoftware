#ifndef FUN4ALLRAW_SINGLEMVTXPOOLINPUT_H
#define FUN4ALLRAW_SINGLEMVTXPOOLINPUT_H

#include "SingleStreamingInput.h"

#include <algorithm>
#include <map>
#include <vector>

class MvtxRawHit;
class Packet;

typedef struct linkId
{
  uint32_t layer = 0xFF;
  uint32_t stave = 0xFF;
  uint32_t gbtid = 0xFF;
} LinkId_t;

class SingleMvtxPoolInput : public SingleStreamingInput
{
 public:
  explicit SingleMvtxPoolInput(const std::string &name);
  ~SingleMvtxPoolInput() override;
  void FillPool(const unsigned int nevents = 1) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  void SetBcoRange(const unsigned int i) { m_BcoRange = i; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) {m_NegativeBco = value;}

  std::set<int> &getFeeIdSet(const uint64_t &bco) { return m_BeamClockFEE[bco]; };

 protected:
  LinkId_t DecodeFeeid(const uint16_t &feeid)
  {
    LinkId_t ret = {};
    ret.layer = (feeid >> 12) & 0x7;
    ret.stave = feeid & 0x1F;
    ret.gbtid = (feeid >> 8) & 0x3;
    return ret;
  }

 private:
  Packet **plist{nullptr};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco {0};

  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<MvtxRawHit *>> m_MvtxRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::map<int, uint64_t> m_FeeStrobeMap;
  std::set<uint64_t> m_BclkStack;
  std::set<uint64_t> gtmL1BcoSet;  // GTM L1 BCO
};

#endif
