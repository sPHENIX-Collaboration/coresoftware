#ifndef GLOBALQA_GLOBALQA_H
#define GLOBALQA_GLOBALQA_H

#include <fun4all/SubsysReco.h>

#include <limits>  // for numeric_limits
#include <string>
#include <vector>

// Forward declarations
class CDBTTree;
class PHCompositeNode;
class TH1;
class TH2;
class TProfile2D;

class GlobalQA : public SubsysReco
{
 public:
  //! constructor
  GlobalQA(const std::string &name =
               "GlobalQA");  // const std::string &filename = "testQA.root");
  // //int nevents = 100);

  //! destructor
  ~GlobalQA() override;

  //! full initialization
  int Init(PHCompositeNode *) override;

  //! event processing method
  int process_event(PHCompositeNode *) override;

  int process_g4hits(PHCompositeNode *);
  int process_g4cells(PHCompositeNode *);
  int process_towers(PHCompositeNode *);
  int process_clusters(PHCompositeNode *);

  void Detector(const std::string &name) { detector = name; }
  void set_timing_cut_width(const int &t) { _range = t; }

  void set_debug(bool debug) { m_debug = debug; }
  TH2 *LogYHist2D(const std::string &name, const std::string &title, int,
                  double, double, int, double, double);

 private:
  int evtcount{0};
  int Getpeaktime(TH1 *h);
  void createHistos();

  // MBD histos
  TH1 *h_GlobalQA_mbd_zvtx{nullptr};
  TH1 *h_GlobalQA_mbd_zvtx_wide{nullptr};
  TH1 *h_GlobalQA_calc_zvtx{nullptr};
  TH1 *h_GlobalQA_calc_zvtx_wide{nullptr};
  TH1 *h_GlobalQA_mbd_charge_s{nullptr};
  TH1 *h_GlobalQA_mbd_charge_n{nullptr};
  TH1 *h_GlobalQA_mbd_charge_sum{nullptr};
  TH2 *h2_GlobalQA_mbd_charge_NS_correlation{nullptr};
  TH1 *h_GlobalQA_mbd_nhit_s{nullptr};
  TH1 *h_GlobalQA_mbd_nhit_n{nullptr};
  TH2 *h2_GlobalQA_mbd_nhits_NS_correlation{nullptr};
  TH1 *h_GlobalQA_mbd_zvtxq{nullptr};

  // sEPD
  TH1 *h_GlobalQA_sEPD_tile[744] = {nullptr};
  TH1 *h_GlobalQA_sEPD_adcsum_s{nullptr};
  TH1 *h_GlobalQA_sEPD_adcsum_n{nullptr};
  TH2 *h2_GlobalQA_sEPD_adcsum_ns{nullptr};
  TH2 *h2_GlobalQA_sEPD_ADC_channel_north{nullptr};
  TH2 *h2_GlobalQA_sEPD_ADC_channel_south{nullptr};
  TProfile2D *h2Profile_GlobalQA_sEPD_tiles_north{nullptr};
  TProfile2D *h2Profile_GlobalQA_sEPD_tiles_south{nullptr};

  // ZDC histos
  TH1 *h_GlobalQA_zdc_zvtx{nullptr};
  TH1 *h_GlobalQA_zdc_zvtx_wide{nullptr};
  TH1 *h_GlobalQA_zdc_energy_s{nullptr};
  TH1 *h_GlobalQA_zdc_energy_n{nullptr};
  TH1 *h_GlobalQA_triggerVec{nullptr};
  float zdc_zvtx{std::numeric_limits<float>::quiet_NaN()};

  int _eventcounter{0};
  int _range{1};

  bool m_debug{false};

  std::string detector;
  std::string m_outputFileName;
  std::string OutputFileName;

  // MBD trigger bit definitions (matching online monitoring)
  static constexpr uint64_t mbdns = (0x1UL << 10) | (0x1UL << 11);  // MBD NS2 and NS1
  static constexpr uint64_t mbdnsvtx10 = (0x1UL << 12) | (0x1UL << 15);  // MBD NS2/NS1 with vtx < 10cm
  static constexpr uint64_t mbdnsvtx30 = (0x1UL << 13);  // MBD NS2 with vtx < 30cm
  static constexpr uint64_t mbdnsvtx150 = (0x1UL << 14);  // MBD NS2 with vtx < 150cm
  static constexpr uint64_t mbdtrig = mbdns | mbdnsvtx10 | mbdnsvtx30 | mbdnsvtx150;  // Combined MBD triggers
};

#endif
