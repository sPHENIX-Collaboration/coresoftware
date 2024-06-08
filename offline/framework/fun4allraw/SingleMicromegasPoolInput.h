#ifndef FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H
#define FUN4ALLRAW_SINGLEMICROMEGASPOOLINPUT_H

#include "SingleStreamingInput.h"
#include "MicromegasBcoMatchingInformation.h"

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

  /// enable evaluation
  void SetDoEvaluation(bool value)
  { m_do_evaluation = value; }

  /// output file name for evaluation histograms
  void SetEvaluationOutputFilename(const std::string& outputfile)
  { m_evaluation_filename = outputfile; }

  //! save some statistics for BCO statistics
  void FillBcoStatistics( uint64_t /*gtm_bco*/);

 private:
  std::array<Packet*,10> plist{};
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

  // keep track of total number of waveforms
  uint64_t m_waveform_count_total = 0;

  // keep track of dropped waveforms
  uint64_t m_waveform_count_dropped = 0;

  bool m_do_evaluation = false;
  std::string m_evaluation_filename = "SingleMicromegasPoolInput.root";
  std::unique_ptr<TFile> m_evaluation_file;

  //!@name gtm bco statistics histogram
  //@{
  TH1* m_npacket_bco_hist = nullptr;
  TH1* m_nwaveform_bco_hist = nullptr;
  //@}

};

#endif
