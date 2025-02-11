// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QA_INTT_INTTRAWHITQA_H
#define QA_INTT_INTTRAWHITQA_H

// Fun4All libraries
#include <fun4all/SubsysReco.h>

// std libraries
#include <string>
#include <vector>

class InttRawHit;
class InttRawHitContainer;
class PHCompositeNode;
class TH1;
class TH2;
class TH3;

class InttRawHitQA : public SubsysReco
{
 public:
  InttRawHitQA(const std::string& name = "InttRawHitQA");

  ~InttRawHitQA() override = default;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode* topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode* topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode* topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode* topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode* topNode) override;

  /// Reset
  int Reset(PHCompositeNode* /*topNode*/) override;

 private:
  void createHistos();
  std::string getHistoPrefix() const;
  ///////////////////////////////////////////
  // general variables
  ///////////////////////////////////////////
  static const int kFelix_num_ = 8;     // the number of our FELIX server
  static const int kFee_num_ = 14;      // the number of half-ladders in a single FELIX server
  static const int kChip_num_ = 26;     // the number of chip in a half-ladder
  static const int kChan_num_ = 128;    // the number of channel in a single chip
  static const int kFirst_pid_ = 3001;  // the first pid (packet ID), which means intt0

    std::vector<InttRawHitContainer*> m_rawhit_containers;
  int previous_event_counter_ = -1;
  int last_event_counter_ = 0;
  int event_counter_by_myself_ = 0;  // because the event counter is not reliable, I count it by myself for histogram normalization

  bool is_first_event_ = true;

  ///////////////////////////////////////////
  // objects to be output
  ///////////////////////////////////////////

  // mother 3D hist
  TH3* hist_fee_chip_chan_[kFelix_num_]{nullptr};  // ch vs chip vs ladder vs felix
  // TH3* hist_fee_chip_chan_woclonehit_[ kFelix_num_ ]; // ch vs chip vs ladder vs felix ; without clonehit
  TH3* hist_fee_bco_full_event_counter_[kFelix_num_]{nullptr};       // event counter vs bco full vs ladder vs felix
  TH3* hist_fee_bco_full_event_counter_diff_[kFelix_num_]{nullptr};  // difference of event counter vs difference of bco full vs ladder vs felix, difference means ( val - Min( val(felix=Min(felix) ) ) )

  // 2D hists
  TH2* hist_hitmap_[kFelix_num_][kFee_num_]{{nullptr}};

  // a simple 1D hists
  TH1* hist_nhit_{nullptr};        // the number of INTTRAWHIT
  TH1* hist_pid_{nullptr};         // the number of hits for each FELIX server
  TH1* hist_nhit_south_{nullptr};  // the number of INTTRAWHIT
  TH1* hist_nhit_north_{nullptr};  // the number of INTTRAWHIT

  // TH1* hist_fee_;
  // TH1* hist_chip_;
  // TH1* hist_chan_;
  TH1* hist_adc_{nullptr};
  TH1* hist_bco_{nullptr};       // FPHX BCO
  TH1* hist_bco_full_{nullptr};  // BCO full

  // felix vs event counter
  TH1* hist_event_counter_[kFelix_num_]{nullptr};
  TH1* hist_event_counter_diff_[kFelix_num_]{nullptr};

  ///////////////////////////////////////////
  // functions
  ///////////////////////////////////////////
  virtual std::vector<InttRawHit*> GetHits(InttRawHitContainer* container);
};

#endif  // INTTRAWHITQA_H
