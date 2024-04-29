
// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOTOWERSLOPE_LITECALOEVAL_H
#define CALOTOWERSLOPE_LITECALOEVAL_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TH1;
class TH2;
class TH3;
class TGraph;
class TNtuple;
class TF1;

double LCE_fitf(double *f, double *p);

TGraph *LCE_grff{nullptr};

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

  LiteCaloEval(const std::string &name = "LiteCaloEval", const std::string &caloNm = "CEMC", const std::string &fnm = "outJF");

  // to distinguish when we want to implement input decal (for simulations work)
  void set_mode(int modeset)
  {
    mode = modeset;
  }

  void set_UseTowerInfo(int setTowerInfo)
  {
    m_UseTowerInfo = setTowerInfo;
  }

  virtual ~LiteCaloEval() {}

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place

      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void CaloType(const Calo i)
  {
    calotype = i;
  }

  TFile *f_temp{nullptr};

  void Get_Histos(const std::string &infile, const std::string &fun4all_file = "");

  void FitRelativeShifts(LiteCaloEval *ref_lce, int modeFitShifts);

  /// Setters
  void setFitMax(float fitMax) { fitmax = fitMax; }
  void setFitMin(float fitMin) { fitmin = fitMin; }
  void set_spectra_binWidth(double binWidth) { binwidth = binWidth; }

  bool chk_isChimney(int, int);

  /// Getters
  float getFitMax() { return fitmax; }
  float getFitMin() { return fitmin; }

  void setInputTowerNodeName(const std::string &inpNodenm)
  {
    _inputnodename = inpNodenm;
  }

 private:
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

  Calo calotype{NONE};
  int _ievent{0};

  float fitmin{0.};
  float fitmax{0.};

  double binwidth{0.001};

  std::string _caloname;
  std::string _filename;
  std::string _inputnodename;

  int mode = 0;

  // flag for using tower info
  int m_UseTowerInfo{0};
};

#endif  // LITECALOEVAL_H
