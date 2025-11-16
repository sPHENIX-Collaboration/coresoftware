#ifndef QA_INTT_INTTMIXUPQA_H
#define QA_INTT_INTTMIXUPQA_H

// Fun4All headers
#include <fun4all/SubsysReco.h>

#include <fstream>
#include <map>
#include <set>

class PHCompositeNode;
class TFile;
class TGraph;
class TH1;
class TH2;

class InttMixupQA : public SubsysReco
{
 public:
  explicit InttMixupQA(const std::string &name = "InttMixupQA", const int run_num = 0, const int felix_num = 0);

  ~InttMixupQA() override = default;

  int Init(PHCompositeNode *) override;

  int InitRun(PHCompositeNode *) override;

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  int SetOutputDir(std::string const &dir);

  void SetBcoPeakFileDir(std::string const &path) { bcopeak_dir_ = path; };

  void SetHotChanFileDir(std::string const &path) { hotchan_dir_ = path; };

  std::string getHistoPrefix() { return "InttMixupQA"; }
  // int SetHistBin(std::string type);
 private:
  // general variables
  int run_num_{0};
  int felix_num_{0};
  static const int kFelix_num_{8};     // the number of our FELIX server
  static const int kFee_num_{14};      // the number of half-ladders in a single FELIX server
  static const int kChip_num_{26};     // the number of chip in a half-ladder
  static const int kChan_num_{128};    // the number of channel in a single chip
  static const int kFirst_pid_{3001};  // the first pid (packet ID), which means intt0
  static const int divimul{10};

  // variables for the output
  std::string output_dir_{"./"};
  std::string output_basename_{"InttMixupEventQA_run"};
  std::string output_root_{"mixup.root"};
  std::string output_pdf_{"mixup.pdf"};
  std::string output_txt_{"bcopeak.pdf"};
  TFile *tf_output_{nullptr};

  // variables for get bco peak
  std::string bcopeak_dir_{"./"};
  // std::string bcopeak_file;
  std::string bcopeak_file[kFelix_num_];
  bool force_suffix_{false};
  bool is_official_{false};
  // bool is_official_=true;
  std::string suffix_;

  std::string GetFileSuffix()
  {
    if (force_suffix_ == false)
    {
      if (is_official_ == true)
        return "_official";
      else
        return "_special";
    }
    return suffix_;
  };

  // variables for hot channel cut
  std::string hotchan_dir_{"./"};
  std::string hotchan_file;

  std::set<int> bcopar_[kFelix_num_];
  std::set<int> otbcopar_[kFelix_num_];
  std::map<int, int> hotmap;

  int ievent_{0};
  int n{kFelix_num_};

  int prev_bcofull{0};
  uint64_t long_prev_bcofull{0};
  int pre_allhit[kFelix_num_]{};

  int NmixupEv[kFelix_num_]{};
  double mixupfraction[kFelix_num_]{};
  double mixupfraction_sum[kFelix_num_]{};

  int Nmixup_sum[kFelix_num_]{};
  int pre_allhit_sum[kFelix_num_]{};
  double Nmixup_ave[kFelix_num_]{};
  double pre_allhit_ave[kFelix_num_]{};

  // double thisclonefraction[kFelix_num_];
  double copyfraction[kFelix_num_]{};
  double copyfraction_sum[kFelix_num_]{};
  int NmixcopyEv[kFelix_num_]{};

  double mixupfraction_ave[kFelix_num_]{};
  double Mixevent[kFelix_num_]{};
  double err[kFelix_num_]{};
  double copyfraction_ave[kFelix_num_]{};
  double Mixcloevent[kFelix_num_]{};
  double thisclone_ave[kFelix_num_]{};

  uint64_t first_bcofull{0};

  int DEFAULT_BCO_VALUE[kFelix_num_]{};

  TGraph *g_evfraction{nullptr};
  TGraph *g_cloevfraction{nullptr};
  TGraph *g_hitfraction{nullptr};
  TGraph *g_copyfraction{nullptr};

  std::map<int, int> premap_hit;

  TH1 *h_allmulti_[kFelix_num_]{};
  TH1 *h_allclone_[kFelix_num_]{};

  TH1 *h_mixupmulti[kFelix_num_][divimul]{};
  TH1 *h_divmul[kFelix_num_][divimul]{};

  // plot these
  TH2 *h_vsprefull_bco[kFelix_num_]{};
  TH2 *h_vsfull_bco[kFelix_num_]{};
  TH1 *h_prefull_bco[kFelix_num_]{};
  TH1 *h_full_bco[kFelix_num_]{};

  TH1 *h_prefull_bco_all[kFelix_num_]{};
  TH2 *h_vsprefull_bco_all[kFelix_num_]{};
  TH1 *h_interval{nullptr};
  TH1 *h_mixinterval{nullptr};
  TH1 *h_divinter{nullptr};
  TH1 *h_bcofull_7{nullptr};
  TH1 *h_bco{nullptr};

  TH1 *h_mixup[kFelix_num_]{};
  TH2 *h_prevsNmix[kFelix_num_]{};
  TH1 *h_nocopyhit[kFelix_num_]{};
  TH1 *h_copyhit[kFelix_num_]{};
  TH1 *h_mixcopy[kFelix_num_]{};
  TH2 *h_mixvscopy[kFelix_num_]{};

  TFile *tf_bcopeak_[kFelix_num_]{};
  TH1 *hbco[kFelix_num_]{};
  TH1 *hbco_sub[kFelix_num_]{};
  TH2 *h2_bco_felix[kFelix_num_]{};
  TH2 *h2_bco_felix_sub[kFelix_num_]{};
  TH1 *hbcohist[kFelix_num_]{};
  TH1 *hbcohist2[kFelix_num_]{};
  // TFile *tf_hotchan_;

  // Mixup fraction hist
  TH2 *h_hitfra[kFelix_num_]{};
  TH2 *h_bghit[kFelix_num_]{};
  TH1 *h_NmixEv{nullptr};
  TH1 *h_AllEv{nullptr};

  std::ofstream f_felixpeak;
  std::ofstream ffraction;
  std::ofstream fNhit[kFelix_num_];
  std::ofstream f_hotchan;
  TFile *fgraph{nullptr};

  // for Histgram parameter
  // int bin {1000};
  // int bin {400};
  int bin{8000};
  int bin3{100};
  uint64_t bit{0xFFFFFFFFFF};
  int64_t bin2{10000};
  // std::string p;
  // std::string Au;

  void DrawHists();

  // void Mixupfraction();

  void GetBcopeak();

  void Readpeak();

  void GetBcomyself();

  void Hotchancut();
};
#endif
