#ifndef RAWBCOLUMI_SINGLEGL1POOLINPUTV2_H
#define RAWBCOLUMI_SINGLEGL1POOLINPUTV2_H

#include "SingleStreamingInputv2.h"

#include <cstdint>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class Gl1Packet;
class PHCompositeNode;

class SingleGl1PoolInputv2 : public SingleStreamingInputv2
{
 public:
  explicit SingleGl1PoolInputv2(const std::string &name);
  ~SingleGl1PoolInputv2() override;
  void FillPool(const unsigned int) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;
  void SetBcoRange(const unsigned int i) { m_BcoRange = i; }
  //  void ConfigureStreamingInputManager() override;
  void SetNegativeWindow(const unsigned int value) { m_negative_bco_window = value; }
  void SetPositiveWindow(const unsigned int value) { m_positive_bco_window = value; }
  void SetTotalEvent(const int value) { m_total_event = value; }

 private:
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};

  //! map bco to packet
  std::map<unsigned int, uint64_t> m_packet_bco;

  std::map<uint64_t, std::vector<Gl1Packet *>> m_Gl1RawHitMap;
  std::map<uint64_t, std::pair<uint64_t, uint64_t>> m_BCOWindows;
  std::map<uint64_t, int> m_BCOBunchNumber;
  std::set<uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;

  unsigned int m_negative_bco_window = 20;
  unsigned int m_positive_bco_window = 325;
  bool m_alldone_flag = {false};
  bool m_lastevent_flag = {false};
  int m_total_event = std::numeric_limits<int>::max();
};

#endif
