#ifndef FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_V2_H
#define FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_V2_H

#include "MicromegasBcoMatchingInformation_v2.h"
#include "SingleStreamingInput.h"

#include <phool/PHTimer.h>

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
class TTree;
class TH1;

class SingleMicromegasPoolInput_v2 : public SingleStreamingInput
{
 public:

  /// constructor
  explicit SingleMicromegasPoolInput_v2(const std::string &name = "SingleMicromegasPoolInput_v2");

  /// destructor
  ~SingleMicromegasPoolInput_v2() override;

  /// pool filling
  void FillPool(const uint64_t /*target_bco*/) override;

  /// cleanup
  void CleanupUsedPackets(const uint64_t bclk) override
  {
    CleanupUsedPackets(bclk, false);
  }

  /// specialized verion of cleaning up packets, with an extra flag about wheter the cleanup hits are dropped or not
  void CleanupUsedPackets(const uint64_t /* bclk */, bool /*dropped */) override;

  /// current event cleaning
  void ClearCurrentEvent() override;

  /// print
  void Print(const std::string &what = "ALL") const override;

  ///
  void CreateDSTNode(PHCompositeNode *topNode) override;

  void SetRecoverTruncatedWaveforms( bool value ) { m_recover_truncated_waveforms = value; }

  void SetBcoRange(const unsigned int value) { m_BcoRange = value; }

  void ConfigureStreamingInputManager() override;

  void SetNegativeBco(const unsigned int value) { m_NegativeBco = value; }

  /// set the offset between FELIX tagger BCO (internal) and GL1 (external) BCO
  /**
   * explicitly taggerBCO = GL1 BCO + offset
   * this is somewhat redundant with m_Negative BCO, unfortunately, but
   * 1/ m_NegativeBco is also used upstream by Fun4AllStreamingInputManager
   * 2/ m_NegativeBCO is unsigned int
  */
  void SetTaggerBcoOffset( const int value ) { m_TaggerBcoOffset = value; }

  /// define minimum pool size in terms of how many BCO are stored
  /** deprecated */
  void SetBcoPoolSize(const unsigned int /*value*/)
  { std::cout << "SingleMicromegasPoolInput_v2::SetBcoPoolSize is deprecated" << std::endl; }

  /// save some statistics for BCO QA
  void FillBcoQA(uint64_t /*gtm_bco*/) override;

  // write the initial histograms for QA manager
  void createQAHistos() override;

  /// do evalutation
  void set_do_evaluation(bool value) { m_do_evaluation = value; }

  /// output file name for evaluation histograms
  void set_evaluation_outputfile(const std::string &outputfile) { m_evaluation_filename = outputfile; }

  private:

  /// true if more data is to be processed for collecting that of a given bco
  bool is_more_data_required(const uint64_t /*target_bco*/) const;

  ///@name decoding constants
  //@{
  /// max number of FEE per OBDC
  static constexpr uint16_t MAX_FEECOUNT = 26;

  /// max number of channels per FEE
  static constexpr uint16_t MAX_FEECHANNELCOUNT = 256;

  // Length for the 256-bit wide Round Robin Multiplexer for the data stream
  static constexpr size_t DAM_DMA_WORD_LENGTH = 16;
  //@}

  /// DMA word structure
  struct dma_word
  {
    uint16_t dma_header;
    uint16_t data[DAM_DMA_WORD_LENGTH - 1];
  };

  void process_packet(Packet *);
  void decode_gtm_data(int /*packet_id*/, const dma_word &);
  void process_fee_data(int /*packet_id*/, unsigned int /*fee_id*/);

  /// fill evaluation tree
  void fill_evaluation_tree( const uint64_t /*target_bco*/ );

  /// recover truncated waveforms for a given gtm bco
  void recover_truncated_waveforms( const uint64_t /*target_bco*/ );

  /// true to recover waveform truncated due to overlapping timeframes
  bool m_recover_truncated_waveforms{true};

  /// fee data buffer
  std::array<std::deque<uint16_t>, MAX_FEECOUNT> m_feeData{};

  /// list of packets from data stream
  std::array<Packet *, 10> plist{};

  /// keep track of number of non data events
  unsigned int m_NumSpecialEvents{0};

  /// bco adjustment for matching across subsystems
  unsigned int m_BcoRange{0};

  /// bco adjustment for matching across subsystems
  unsigned int m_NegativeBco{0};

  /// offset between FELIX tagger BCO (internal) and GL1 (external) BCO
  int m_TaggerBcoOffset{3};

  /// store list of packets that have data for a given beam clock
  /**
   * all packets in taggers are stored,
   * disregarding whether there is data associated to it or not
   * this allows to keep track of dropped data, also in zero-suppression mode
   */
  std::map<uint64_t, std::set<int>> m_BeamClockPacket;

  /// store list of FEE that have data for a given beam clock
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;

  /// list of raw hits
  using rawhit_list_t = std::vector<MicromegasRawHit*>;

  /// maps list of raw hits on GTM BCO values
  using rawhit_map_t = std::map<uint64_t, rawhit_list_t>;

  /// store list of raw hits matching a given GTM bco on a per FEE basis
  std::array<rawhit_map_t,MAX_FEECOUNT> m_MicromegasRawHitMap{};

  /// map bco_information_t to packet id
  using bco_matching_information_map_t = std::map<unsigned int, MicromegasBcoMatchingInformation_v2>;
  bco_matching_information_map_t m_bco_matching_information_map{};

  /// map packet to FEE ID
  /* it is filled on the fly. It allows to quickly retrieve BCO matching information from FEE index */
  std::array<unsigned int,MAX_FEECOUNT> m_fee_packet{};

  class counter_t
  {
   public:
    /// total count
    uint64_t total {0};

    /// drop count due to unmatched bco
    uint64_t dropped_bco {0};

    /// drop count due to pools
    uint64_t dropped_pool {0};

    /// dropped fraction (bco)
    double dropped_fraction_bco() const { return double(dropped_bco) / total; }

    /// dropped fraction (pool)
    double dropped_fraction_pool() const { return double(dropped_pool) / total; }
  };

  // keep track of waveform statistics per fee
  std::map<int, counter_t> m_fee_waveform_counters{};

  // keep track of waveform statistics per packet
  std::map<int, counter_t> m_waveform_counters{};

  // keep track of heartbeat statistics per fee
  std::map<int, counter_t> m_fee_heartbeat_counters{};

  // keep track of heartbeat statistics per packet
  std::map<int, counter_t> m_heartbeat_counters{};

  // timer
  PHTimer m_timer{"SingleMicromegasPoolInput_v2"};

  ///@name QA histograms
  //@{

  /// keeps track of how often a given (or all) packets are found for a given BCO
  TH1 *h_packet_stat{nullptr};

  /// keep track of how many heartbeats are found per FEE sampa
  TH1 *h_heartbeat_stat{nullptr};

  /// keeps track of how many packets are found for a given BCO
  TH1 *h_packet{nullptr};

  /// keeps track of how many waveforms are found for a given BCO
  TH1 *h_waveform{nullptr};

  /// total number of waveforms per packet
  TH1 *h_waveform_count_total{nullptr};

  /// total number of dropped waveforms per packet due to bco mismatch
  /*! waveforms are dropped when their FEE-BCO cannot be associated to any global BCO */
  TH1 *h_waveform_count_dropped_bco{nullptr};

  /// total number of dropped waveforms per packet due to fun4all pool mismatch
  TH1 *h_waveform_count_dropped_pool{nullptr};

  /// total number of waveforms per packet
  TH1 *h_fee_waveform_count_total{nullptr};

  /// total number of dropped waveforms per fee due to bco mismatch
  /*! waveforms are dropped when their FEE-BCO cannot be associated to any global BCO */
  TH1 *h_fee_waveform_count_dropped_bco{nullptr};

  /// total number of dropped waveforms per fee due to fun4all pool mismatch
  TH1 *h_fee_waveform_count_dropped_pool{nullptr};

  //@}

  ///@name evaluation
  //@{

  /// evaluation
  bool m_do_evaluation = false;

  /// evaluation output filename
  std::string m_evaluation_filename = "SingleMicromegasPoolInput.root";
  std::unique_ptr<TFile> m_evaluation_file;

  /**
   * waveform is similar to sample except that there is only one of which per waveform,
   * and that it stores the max adc and corresponding sample_id
   */
  class Waveform
  {
   public:
    /// packet
    unsigned int packet_id {0};

    /// fee
    unsigned short fee_id {0};

    /// channel id
    unsigned short channel {0};

    /// ll1 bco
    uint64_t gtm_bco_first {0};

    /// bco
    uint64_t gtm_bco_tagger {0};

    /// bco
    uint64_t gtm_bco_gl1 {0};

    /// fee bco
    unsigned int fee_bco_first {0};

    /// fee bco
    unsigned int fee_bco {0};

    /// fee bco predicted (from gtm tagger)
    unsigned int fee_bco_predicted_tagger {0};

    /// fee bco predicted (from gtm gl1)
    unsigned int fee_bco_predicted_gl1 {0};

  };

  Waveform m_waveform;

  /// tree
  TTree *m_evaluation_tree {nullptr};

  //*}
};

#endif
