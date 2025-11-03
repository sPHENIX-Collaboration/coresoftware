// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HCALCOSMICS_H
#define HCALCOSMICS_H

#include <fun4all/SubsysReco.h>

#include <string>

// Forward declarations
class PHCompositeNode;
class TFile;
class TF1;
class TH1;
class TH2;

class HCalCosmics : public SubsysReco
{
 public:
  //! constructor
  HCalCosmics(const std::string &, const std::string &);
  //! destructor
  ~HCalCosmics() override = default;

  //! Processing
  int Init(PHCompositeNode *) override;
  int process_event(PHCompositeNode *) override;
  int process_towers(PHCompositeNode *);
  int End(PHCompositeNode *) override;

  void set_tower_threshold(float fac) { tower_threshold = fac; }
  void set_vert_threshold(float fac) { vert_threshold = fac; }
  void set_veto_threshold(float fac) { veto_threshold = fac; }

  void set_adc_tower_threshold(float fac) { adc_tower_threshold = fac; }
  void set_adc_vert_threshold(float fac) { adc_vert_threshold = fac; }
  void set_adc_veto_threshold(float fac) { adc_veto_threshold = fac; }

  void HistBinWidth(double fac) { bin_width = fac; }
  void Detector(const std::string &name) { detector = name; }
  void TowerPrefix(const std::string &name) { prefix = name; }

  void RawHistBinWidth(double fac) { rawbin_width = fac; }
  void RawDetector(const std::string &name) { rawdetector = name; }
  void RawTowerPrefix(const std::string &name) { rawprefix = name; }

  static double gamma_function(const double *x, const double *par);
  void fitChannels(const std::string &infile, const std::string &outfilename2);
  static TF1 *fitHist(TH1 *);

 private:
  // HCal geometry
  static const int n_etabin{24};
  static const int n_phibin{64};

  TFile *outfile{nullptr};
  TH1 *h_channel_hist[n_etabin][n_phibin]{{nullptr}};   // distribution of tower energy
  TH1 *h_adc_hist[n_etabin][n_phibin]{{nullptr}};       // distribuiton of raw tower ADC with topo cuts in ADC
  TH1 *h_adc_ecut_hist[n_etabin][n_phibin]{{nullptr}};  // distribution of raw tower ADC with topo cuts in energy
  TH1 *h_gain_hist[n_etabin][n_phibin]{{nullptr}};      // distribution of gain per tower
  TH2 *h_waveformchi2{nullptr};
  TH2 *h_waveformchi2_aftercut{nullptr};
  TH1 *h_towertime_hist[n_etabin][n_phibin]{{nullptr}};  // distribution of tower time of muons
  TH2 *h_time_energy{nullptr};
  TH1 *h_mip{nullptr};
  TH1 *h_adc{nullptr};
  TH1 *h_adc_ecut{nullptr};
  TH1 *h_gain{nullptr};
  TH1 *h_event{nullptr};

  // Cut threshold in GeV
  float tower_threshold{0.2498};  // 500 ADC  iHCal: 0.2498  oHCal: 1.665
  float vert_threshold{0.2498};   // 500 ADC  iHCal: 0.2498  oHCal: 1.665
  float veto_threshold{0.17486};  // 350 ADC  iHCal: 0.17486  oHCal: 1.1655
  // Cut threshold in ADC
  float adc_tower_threshold{500};
  float adc_vert_threshold{500};
  float adc_veto_threshold{350};

  int event{0};

  float bin_width{0.01};  // 20 ADC  iHCal: 0.01  oHCal: 0.05
  float rawbin_width{1};

  float m_peak[n_etabin][n_phibin]{};
  float m_adc[n_etabin][n_phibin]{};
  float m_chi2[n_etabin][n_phibin]{};
  float m_time[n_etabin][n_phibin]{};

  //  bool debug {false};

  std::string prefix{"TOWERINFO_CALIB_"};  // for calibrated tower energy
  std::string detector{"HCALOUT"};

  std::string rawprefix{"TOWERS_"};  // for raw tower ADC
  std::string rawdetector{"HCALOUT"};

  std::string outfilename;
};

#endif  // HCALCOSMICS_H
