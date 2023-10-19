#ifndef FUN4ALLRAW_SINGLEMVTXINPUT_H
#define FUN4ALLRAW_SINGLEMVTXINPUT_H

#include "SingleStreamingInput.h"

//#include <array>
//#include <list>
#include <map>
#include <algorithm>
//#include <set>
//#include <string>
#include <vector>
//
//class Eventiterator;
//class Fun4AllEvtInputPoolManager;
class MvtxRawHit;
class Packet;

typedef struct linkId
{
  uint32_t layer = 0xFF;
  uint32_t stave = 0xFF;
  uint32_t gbtid = 0xFF;
} LinkId_t;

//struct MvtxL1Trg
//{
//  uint64_t bco : 40; // 40 bits gl1/gtm bco
//  uint16_t bc  : 12; // 12 bits mvtx internal bc counter
//
//  MvtxL1Trg() = default;
//  ~MvtxL1Trg() = default;
//
//  MvtxL1Trg(uint64_t _bco, uint16_t _bc) : bco(_bco), bc(_bc) {};
//
//  bool operator<(const MvtxL1Trg& other) const
//  {
//    return (bco == other.bco) ? (bc < other.bc) : (bco < other.bco);
//  }
//
//  bool operator==(const MvtxL1Trg& other) const
//  {
//    return (bco == other.bco) ? (bc == other.bc) : false;
//  }
//
//};

class SingleMvtxInput : public SingleStreamingInput
{
 public:
  explicit SingleMvtxInput(const std::string &name);
  ~SingleMvtxInput() override;
  void FillPool(const unsigned int nevents = 1) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

  std::set<uint64_t>& getGtmL1BcoSet() {return m_GtmL1BcoSet; };
  std::set<int>& getFeeIdSet(const uint64_t& bco) { return m_BeamClockFEE[bco]; };
  void clearGtmL1BcoSet();

 protected:
  LinkId_t DecodeFeeid(const uint16_t& feeid)
  {
    LinkId_t ret = {};
    ret.layer = (feeid >> 12) & 0x7;
    ret.stave = feeid & 0x1F;
    ret.gbtid = (feeid >> 8) & 0x3;
    return ret;
  }

 private:
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<MvtxRawHit *>> m_MvtxRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::map<int, uint64_t> m_FeeStrobeMap;
  std::set<uint64_t> m_BclkStack;
  std::set<uint64_t> m_GtmL1BcoSet;
};

#endif
