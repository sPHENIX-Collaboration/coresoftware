#ifndef FUN4ALLRAW_SINGLEMVTXINPUT_H
#define FUN4ALLRAW_SINGLEMVTXINPUT_H

#include "SingleStreamingInput.h"

//#include <array>
//#include <list>
#include <map>
//#include <set>
//#include <string>
//#include <vector>
//
//class Eventiterator;
//class Fun4AllEvtInputPoolManager;
class MvtxRawHit;
class Packet;

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

 private:
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;
//  std::array<uint64_t, 14> m_PreviousClock{};
//  std::array<uint64_t, 14> m_Rollover{};
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<MvtxRawHit *>> m_MvtxRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif
