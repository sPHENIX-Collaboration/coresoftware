#ifndef FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H
#define FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H

#include "SingleStreamingInput.h"
#include "MicromegasBcoMatchingInformation.h"

#include <array>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class Fun4AllEvtInputPoolManager;
class MicromegasRawHit;
class Packet;

class SingleMicromegasPoolInput : public SingleStreamingInput
{
 public:
  explicit SingleMicromegasPoolInput(const std::string &name = "SingleMicromegasPoolInput");
  ~SingleMicromegasPoolInput() override;
  void FillPool(const unsigned int nevents = 1) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

  void SetBcoRange(const unsigned int value) { m_BcoRange = value; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }

 private:
  std::array<Packet*,10> plist{};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};

  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<MicromegasRawHit *>> m_MicromegasRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;

  //! map bco_information_t to packet id
  using bco_matching_information_map_t = std::map<unsigned int, MicromegasBcoMatchingInformation>;
  bco_matching_information_map_t m_bco_matching_information_map;

  // keep track of total number of waveforms
  uint64_t m_waveform_count_total = 0;

  // keep track of dropped waveforms
  uint64_t m_waveform_count_dropped = 0;

};

#endif
