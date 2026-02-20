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

  void do_skim_EMCal(bool do_skim, uint16_t threshold) {
    b_do_skim_EMCal = do_skim;
    m_EMC_skim_threshold = threshold;
  }

  void do_skim_HCal(bool do_skim, uint16_t threshold) {
    b_do_skim_HCal = do_skim;
    m_HCal_skim_threshold = threshold;
  }

  void do_skim_sEPD(bool do_skim, uint16_t threshold) {
    b_do_skim_sEPD = do_skim;
    m_sEPD_skim_threshold = threshold;
  }

  void do_skim_ZDC(bool do_skim, uint16_t threshold) {
    b_do_skim_ZDC = do_skim;
    m_ZDC_skim_threshold = threshold;
  }

private:
  uint32_t n_eventcounter{0};
  uint32_t n_skimcounter{0};
  uint32_t n_notowernodecounter{0};

  bool b_do_skim_EMCal{false};
  uint16_t m_EMC_skim_threshold{
      192}; // skim if nchannels >= this many not-instrumented(empty/missing
            // pckt) channels in EMCal

  bool b_do_skim_HCal{false};
  uint16_t m_HCal_skim_threshold{
      192}; // skim if nchannels >= this many not-instrumented(empty/missing
            // pckt) channels in HCal

  bool b_do_skim_sEPD{false};
  uint16_t m_sEPD_skim_threshold{
      192}; // skim if nchannels >= this many not-instrumented(empty/missing
            // pckt) channels in sEPD

  bool b_do_skim_ZDC{false};
  uint16_t m_ZDC_skim_threshold{
      192}; // skim if nchannels >= this many not-instrumented(empty/missing
            // pckt) channels in ZDC
};

#endif // CALOSTATUSSKIMMER_H
