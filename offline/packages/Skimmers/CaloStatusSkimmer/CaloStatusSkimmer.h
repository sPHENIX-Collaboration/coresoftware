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

class CaloStatusSkimmer : public SubsysReco {
public:
  CaloStatusSkimmer(const std::string &name = "CaloStatusSkimmer");

  ~CaloStatusSkimmer() override = default;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void do_skim_EMCal( uint16_t threshold) {
    m_EMC_skim_threshold = threshold;
  }

  void do_skim_HCal( uint16_t threshold) {
    m_HCal_skim_threshold = threshold;
  }

  void do_skim_sEPD( uint16_t threshold) {
    m_sEPD_skim_threshold = threshold;
  }

  void do_skim_ZDC( uint16_t threshold) {
    m_ZDC_skim_threshold = threshold;
  }

private:
  uint32_t n_eventcounter{0};
  uint32_t n_skimcounter{0};
  uint32_t n_notowernodecounter{0};

  // If the threshold is set to 0, then the skimming for that subsystem is disabled. If threshold is > 0, then the event is skimmed if nchannels >= threshold not-instrumented (empty/missing packet) channels in that subsystem.
  uint16_t m_EMC_skim_threshold{192}; 
  // skim if nchannels >= this many not-instrumented (empty/missing packet) channels in EMCal

  uint16_t m_HCal_skim_threshold{192}; 
  // skim if nchannels >= this many not-instrumented (empty/missing packet) channels in HCal

  uint16_t m_sEPD_skim_threshold{1}; 
  // skim if nchannels >= this many not-instrumented (empty/missing packet) channels in sEPD

  uint16_t m_ZDC_skim_threshold{1}; 
  // skim if nchannels >= this many not-instrumented (empty/missing packet) channels in ZDC
};

#endif // CALOSTATUSSKIMMER_H
