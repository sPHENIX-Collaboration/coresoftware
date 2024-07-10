#ifndef FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H
#define FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H

#include "MicromegasBcoMatchingInformation.h"
#include "SingleStreamingInput.h"

#include <array>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

class Fun4AllEvtInputPoolManager;
class MicromegasRawHit;
class Packet;

class TFile;
class TH1;

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

  //! save some statistics for BCO QA
  void FillBcoQA(uint64_t /*gtm_bco*/);

  // write the initial histograms for QA manager
  void createQAHistos();

 private:
  std::array<Packet *, 10> plist{};
  unsigned int m_NumSpecialEvents{0};
  unsigned int m_BcoRange{0};
  unsigned int m_NegativeBco{0};

  //! store list of packets that have data for a given beam clock
  /**
   * all packets in taggers are stored,
   * disregarding whether there is data associated to it or not
   * this allows to keep track of dropped data, also in zero-suppression mode
   */
  std::map<uint64_t, std::set<int>> m_BeamClockPacket;

  //! store list of FEE that have data for a given beam clock
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;

  //! store list of raw hits matching a given bco
  std::map<uint64_t, std::vector<MicromegasRawHit *>> m_MicromegasRawHitMap;

  //! store current list of BCO on a per fee basis.
  /** only packets for which a given FEE have data are stored */
  std::map<int, uint64_t> m_FEEBclkMap;

  //! store current list of BCO
  /**
   * all packets in taggers are stored,
   * disregarding whether there is data associated to it or not
   * this allows to keep track of dropped data, also in zero-suppression mode
   */
  std::set<uint64_t> m_BclkStack;

  //! map bco_information_t to packet id
  using bco_matching_information_map_t = std::map<unsigned int, MicromegasBcoMatchingInformation>;
  bco_matching_information_map_t m_bco_matching_information_map;

  // keep track of total number of waveforms per packet
  std::map<int,uint64_t> m_waveform_count_total{};

  // keep track of dropped waveforms per packet
  std::map<int,uint64_t> m_waveform_count_dropped{};

  // QA histograms
  TH1 *h_packet_stat{nullptr};
  TH1 *h_packet{nullptr};
  TH1 *h_waveform{nullptr};
};

#endif
