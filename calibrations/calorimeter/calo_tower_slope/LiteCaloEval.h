// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef LITECALOEVAL_H
#define LITECALOEVAL_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TH1;
class TH2;
class TH3;

class LiteCaloEval : public SubsysReco
{
 public:
  enum Calo
  {
    NONE = 0,
    CEMC = 1,
    HCALIN = 2,
    HCALOUT = 3
  };

  LiteCaloEval(const std::string &name = "LiteCaloEval", const std::string &caloNm = "CEMC", const std::string &fnm = "outJF");

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

  void CaloType(const Calo i) { calotype = i; }

 private:
  TFile *cal_output = nullptr;

  TH1 *hcal_out_eta_phi[24][64] = {};
  TH1 *hcalout_eta[24] = {};
  TH2 *hcalout_energy_eta = nullptr;
  TH3 *hcalout_e_eta_phi = {};

  TH1 *hcal_in_eta_phi[24][64] = {};
  TH1 *hcalin_eta[24] = {};
  TH2 *hcalin_energy_eta = nullptr;
  TH3 *hcalin_e_eta_phi = nullptr;

  TH1 *cemc_hist_eta_phi[96][258] = {};
  TH1 *eta_hist[96] = {};
  TH2 *energy_eta_hist = nullptr;
  TH3 *e_eta_phi = nullptr;

  Calo calotype = NONE;
  int _ievent = 0;
  std::string _caloname;
  std::string _filename;
};

#endif  // LITECALOEVAL_H
