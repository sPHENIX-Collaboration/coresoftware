#ifndef FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_V2_H
#define FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_V2_H

#include "MicromegasBcoMatchingInformation_v2.h"
#include "SingleStreamingInput.h"

#include <array>
#include <deque>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

class MicromegasRawHit;
class Packet;

class TFile;
class TH1;

class SingleMicromegasPoolInput_v2 : public SingleStreamingInput
{
 public:
  explicit SingleMicromegasPoolInput_v2(const std::string &name = "SingleMicromegasPoolInput_v2");
  ~SingleMicromegasPoolInput_v2() override;
  void FillPool(const unsigned int nevents = 1) override;

  void CleanupUsedPackets(const uint64_t bclk) override
  { CleanupUsedPackets_with_qa(bclk,false); }

  //! specialized verion of cleaning up packets, with an extra flag about wheter the cleanup hits are dropped or not
  void CleanupUsedPackets_with_qa(const uint64_t bclk, bool /*dropped */);

  void ClearCurrentEvent() override;
  bool GetSomeMoreEvents();
  void Print(const std::string &what = "ALL") const override;
  void CreateDSTNode(PHCompositeNode *topNode) override;

  void SetBcoRange(const unsigned int value) { m_BcoRange = value; }
  void ConfigureStreamingInputManager() override;
  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }

  //! define minimum pool size in terms of how many BCO are stored
  void SetBcoPoolSize(const unsigned int value) { m_BcoPoolSize = value; }

  //! save some statistics for BCO QA
  void FillBcoQA(uint64_t /*gtm_bco*/);

  // write the initial histograms for QA manager
  void createQAHistos();

 private:

  //!@name decoding constants
  //@{
  // Length for the 256-bit wide Round Robin Multiplexer for the data stream
  static constexpr size_t DAM_DMA_WORD_LENGTH = 16;

  static constexpr uint16_t FEE_PACKET_MAGIC_KEY_1 = 0xfe;
  static constexpr uint16_t FEE_PACKET_MAGIC_KEY_2 = 0xed;

  static constexpr uint16_t FEE_MAGIC_KEY = 0xba00;
  static constexpr uint16_t GTM_MAGIC_KEY = 0xbb00;
  static constexpr uint16_t GTM_LVL1_ACCEPT_MAGIC_KEY = 0xbbf0;
  static constexpr uint16_t GTM_ENDAT_MAGIC_KEY = 0xbbf1;
  static constexpr uint16_t GTM_MODEBIT_MAGIC_KEY = 0xbbf2;

  static constexpr uint16_t MAX_FEECOUNT = 26;      // that many FEEs
  static constexpr uint16_t MAX_SAMPA = 8;      // that many FEEs
  static constexpr uint16_t MAX_CHANNELS = MAX_SAMPA * 32;  // that many channels per FEE
                                                //  static const uint16_t  HEADER_LENGTH  = 5;
  static constexpr uint16_t HEADER_LENGTH = 7;
  static constexpr uint16_t MAX_PACKET_LENGTH = 1025;


  static constexpr size_t MAX_SAMPLE=1024;
  static constexpr int ADC_INVALID=65000;

  //@}

  //! DMA word structure
  struct dma_word
  {
    uint16_t dma_header;
    uint16_t data[DAM_DMA_WORD_LENGTH - 1];
  };

  void process_packet(Packet*);
  void decode_gtm_data(int /*packet_id*/, const dma_word&);
  void process_fee_data(int /*packet_id*/, unsigned int /*fee_id*/);

  // fee data buffer
  std::vector<std::deque<uint16_t>> m_feeData{MAX_FEECOUNT};

  // list of packets from data stream
  std::array<Packet *, 10> plist{};

  /// keep track of number of non data events
  unsigned int m_NumSpecialEvents{0};

  /// bco adjustment for matching across subsystems
  unsigned int m_BcoRange{0};

  /// bco adjustment for matching across subsystems
  unsigned int m_NegativeBco{0};

  //! minimum number of BCO required in Micromegas Pools
  unsigned int m_BcoPoolSize{1};

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
  using bco_matching_information_map_t = std::map<unsigned int, MicromegasBcoMatchingInformation_v2>;
  bco_matching_information_map_t m_bco_matching_information_map;

  // keep track of total number of waveforms per fee
  std::map<int,uint64_t> m_fee_waveform_count_total{};

  // keep track of dropped waveforms per fee, due to BCO mismatched
  std::map<int,uint64_t> m_fee_waveform_count_dropped_bco{};

  // keep track of total number of waveforms per packet
  std::map<int,uint64_t> m_waveform_count_total{};

  // keep track of dropped waveforms per packet, due to BCO mismatched
  std::map<int,uint64_t> m_waveform_count_dropped_bco{};

  // keep track of dropped waveforms per packet, due to fun4all pool mismatch
  std::map<int,uint64_t> m_waveform_count_dropped_pool{};

  //!@name QA histograms
  //@{

  //! keeps track of how often a given (or all) packets are found for a given BCO
  TH1 *h_packet_stat{nullptr};

  //! keeps track of how many packets are found for a given BCO
  TH1 *h_packet{nullptr};

  //! keeps track of how many waveforms are found for a given BCO
  TH1 *h_waveform{nullptr};

  //! total number of waveforms per packet
  TH1 *h_waveform_count_total{nullptr};

  //! total number of dropped waveforms per packet due to bco mismatch
  /*! waveforms are dropped when their FEE-BCO cannot be associated to any global BCO */
  TH1 *h_waveform_count_dropped_bco{nullptr};

  //! total number of dropped waveforms per packet due to fun4all pool mismatch
  TH1 *h_waveform_count_dropped_pool{nullptr};

  //@}

};

#endif
