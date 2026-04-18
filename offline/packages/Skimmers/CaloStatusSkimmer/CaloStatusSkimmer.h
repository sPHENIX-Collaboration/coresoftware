// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOSTATUSSKIMMER_H
#define CALOSTATUSSKIMMER_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <limits>
#include <string>
#include <utility>

class PHCompositeNode;
class TH1;

class CaloStatusSkimmer : public SubsysReco {
public:
  CaloStatusSkimmer(const std::string &name = "CaloStatusSkimmer");

  ~CaloStatusSkimmer() override = default;

  int Init(PHCompositeNode* topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void do_skim_EMCal( uint16_t threshold) 
  {
    m_EMC_skim_threshold = threshold;
  }

  void do_skim_HCal( uint16_t threshold) 
  {
    m_HCal_skim_threshold = threshold;
  }

  void do_skim_sEPD( uint16_t threshold) 
  {
    m_sEPD_skim_threshold = threshold;
  }

  void do_skim_ZDC( uint16_t threshold) 
  {
    m_ZDC_skim_threshold = threshold;
  }

  void produce_QA_histograms(bool produce) 
  {
    b_produce_QA_histograms = produce;
  }

private:
  uint32_t n_eventcounter{0};
  uint32_t n_skimcounter{0};
  uint32_t n_notowernodecounter{0};

  bool b_produce_QA_histograms{false};

  // If the threshold is set to 0, then the skimming for that subsystem is disabled. If threshold is > 0, then the event is skimmed if nchannels >= threshold not-instrumented (empty/missing packet) channels in that subsystem.

  uint16_t m_EMC_skim_threshold{193}; 
  // skim if nchannels >= this many not-instrumented (empty/missing packet) channels in EMCal. For the EMCal in particular we want to skim if greater than 1 packet's worth of channels are not-instrumented, which corresponds to 193 channels (since each packet has 192 channels)

  uint16_t m_HCal_skim_threshold{192}; 
  // skim if nchannels >= this many not-instrumented (empty/missing packet) channels in HCal. Corresponds to 1 packet's worth of channels in HCal, which has 192 channels per packet

  uint16_t m_sEPD_skim_threshold{1}; 
  // skim if nchannels >= this many not-instrumented (empty/missing packet) channels in sEPD. 

  uint16_t m_ZDC_skim_threshold{0}; 
  // skim if nchannels > this many not-instrumented (empty/missing packet) channels in ZDC. Some issue in the ZDC right now so skimming is turned off by now by default.

  // Counters for number of events skimmed per subsystem
  uint32_t EMC_skim_count = 0;
  uint32_t HCal_skim_count = 0;
  uint32_t sEPD_skim_count = 0; 
  uint32_t ZDC_skim_count = 0;

  //Per-calo tower counter histograms
  TH1* h_EMC_nTowers_notinstr = nullptr;
  TH1* h_HCal_nTowers_notinstr = nullptr;
  TH1* h_sEPD_nTowers_notinstr = nullptr;
  TH1* h_ZDC_nTowers_notinstr = nullptr;

  //Event counter histograms
  TH1* h_calo_nEvents = nullptr;

};

#endif // CALOSTATUSSKIMMER_H
