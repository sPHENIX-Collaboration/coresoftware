// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWERSLOPE_LITECALOEVAL_H
#define CALOTOWERSLOPE_LITECALOEVAL_H

#include <fun4all/SubsysReco.h>

#include <string>

class TFile;
class TH1;
class TH2;
class TH3;
class TGraph;
class TriggerAnalyzer;

class LiteCaloEval : public SubsysReco
{
 public:
  int m_myminbin = -1;
  int m_mymaxbin = -3;

  enum Calo
  {
    NONE = 0,
    CEMC = 1,
    HCALIN = 2,
    HCALOUT = 3
  };

  LiteCaloEval(const std::string &name = "LiteCaloEval", const std::string &caloname = "CEMC", const std::string &filename = "outJF");

  virtual ~LiteCaloEval() = default;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void CaloType(const Calo i)
  {
    calotype = i;
  }

  /// Setters_____________________________________________

  void setFitMax(float fitMax) { fitmax = fitMax; }

  void setFitMin(float fitMin) { fitmin = fitMin; }

  void set_spectra_binWidth(double binWidth) { binwidth = binWidth; }

  void set_mode(int modeset)  // to distinguish when we want to implement input decal (for simulations work)
  {
    mode = modeset;
  }

  void set_reqMinBias(bool status)
  {
    reqMinBias = status;
    return;
  }

  void set_doQA(bool status = true)
  {
    doQA = status;
  }

  void setInputTowerNodeName(const std::string &inpNodenm)
  {
    _inputnodename = inpNodenm;
  }

  void set_UseTowerInfo(int setTowerInfo)
  {
    m_UseTowerInfo = setTowerInfo;
  }

  /// Getters________________________________________

  void Get_Histos(const std::string &infile, const std::string &outfile = "");

  float getFitMax() { return fitmax; }

  float getFitMin() { return fitmin; }

  float get_spectra_binWidth() { return binwidth; }

  /// Others__________________________________________

  void FitRelativeShifts(LiteCaloEval *ref_lce, int modeFitShifts);

  static float spec_QA(TH1 *h_spec, TH1 *h_ref, bool retFloat);

  static bool spec_QA(TH1 *h_spec, TH1 *h_ref);

  void plot_cemc(const std::string &path);

  void draw_spectra(const char *);

  void fit_info(const char *, const int);

  void fit_info(const char *, const char *, const int);

  static bool chk_isChimney(int, int);

 private:
  TFile *f_temp{nullptr};
  TFile *cal_output{nullptr};

  TH1 *hcal_out_eta_phi[24][64] = {};
  TH1 *hcalout_eta[25] = {};
  TH2 *hcalout_energy_eta{nullptr};
  TH3 *hcalout_e_eta_phi = {};

  TH1 *hcal_in_eta_phi[24][64] = {};
  TH1 *hcalin_eta[25] = {};
  TH2 *hcalin_energy_eta{nullptr};
  TH3 *hcalin_e_eta_phi{nullptr};

  TH1 *cemc_hist_eta_phi[96][258] = {};
  TH1 *eta_hist[97] = {};
  TH2 *energy_eta_hist{nullptr};
  TH3 *e_eta_phi{nullptr};

  TH1 *h_event{nullptr};

  Calo calotype{NONE};
  int _ievent{0};

  float fitmin{0.};
  float fitmax{0.};

  bool doQA{false};

  double binwidth{0.001};

  std::string _caloname;
  std::string _filename;
  std::string _inputnodename{"TOWERINFO"};

  bool reqMinBias{true};

  int mode{0};

  TriggerAnalyzer *trigAna{nullptr};

  // flag for using tower info
  int m_UseTowerInfo{1};
};

#endif  // LITECALOEVAL_H
