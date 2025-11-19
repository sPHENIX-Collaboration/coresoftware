#ifndef QA_INTT_BCOCHECK_H
#define QA_INTT_BCOCHECK_H

// Fun4All headers
#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TH1;

class bcocheck : public SubsysReco
{
 public:
  explicit bcocheck(const std::string &name = "bcocheck", const int run_num = 0, const int felix_num = 0);

  ~bcocheck() override = default;

  int Init(PHCompositeNode *) override;

  int InitRun(PHCompositeNode *) override;

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  // int SetHistBin(std::string type);
 private:
  // methods
  void DrawHists();
  // general variables
  int run_num_{0};
  int felix_num_{0};
  static constexpr int kFelix_num_{8};     // the number of our FELIX server
  static constexpr int kFee_num_{14};      // the number of half-ladders in a single FELIX server
  static constexpr int kChip_num_{26};     // the number of chip in a half-ladder
  static constexpr int kChan_num_{128};    // the number of channel in a single chip
  static constexpr int kFirst_pid_{3001};  // the first pid (packet ID), which means intt0
  static constexpr int divimul{10};

  int ievent_{0};
  TFile *tf_output_[kFelix_num_]{};

  TH1 *h_full_bco[kFelix_num_]{};
};
#endif
