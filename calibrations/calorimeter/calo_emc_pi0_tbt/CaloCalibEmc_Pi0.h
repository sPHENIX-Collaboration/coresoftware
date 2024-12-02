// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOEMCPI0TBT_CALOCALIBEMC_PI0_H
#define CALOEMCPI0TBT_CALOCALIBEMC_PI0_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>

class PHCompositeNode;
class TFile;
class TH1F;
class TH2F;
class TH3F;
class TH1;
class TTree;

class CaloCalibEmc_Pi0 : public SubsysReco
{
 public:
  CaloCalibEmc_Pi0(const std::string &name = "CaloCalibEmc_Pi0", const std::string &fnm = "outJF");

  ~CaloCalibEmc_Pi0() override = default;

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

  void Loop(int nevts, const std::string &filename, TTree *intree = nullptr, const std::string &ifileCorr = "");
  void Loop_for_eta_slices(int nevts, const std::string &filename, TTree *intree = nullptr, const std::string &ifileCorr = "");

  void Fit_Histos_Etas96(const std::string &infilent);
  void Fit_Histos(const std::string &infilent);
  void Fit_Histos_Eta_Phi_Add96(const std::string &infilent);
  void Fit_Histos_Eta_Phi_Add32(const std::string &infilent);

  void set_centrality_nclusters_cut(int n) { m_cent_nclus_cut = n; }

  void Add_32();
  void Add_96();

  void Get_Histos(const std::string &infile, const std::string &outfile);

  void set_UseTowerInfo(const int useMode)
  {  // 0 only old tower, 1 only new (TowerInfo based),
    m_UseTowerInfo = useMode;
  }

  void setInputClusterNodeName(const std::string &inpNodenm)
  {
    _inputnodename = inpNodenm;
  }

  void setInputTowerNodeName(const std::string &inpNodenm)
  {
    _inputtownodename = inpNodenm;
  }

  void set_calibSetMassVal(float insetval)
  {
    _setMassVal = insetval;
  }

 private:
  //  float setMassVal = 0.135;
  float _setMassVal{0.152};
  // currently defaulting to 0.152 to match sim

  int m_ievent{0};
  std::string m_Filename;
  TFile *cal_output{nullptr};
  std::string _caloname{"CEMC"};
  std::string _inputnodename;
  std::string _inputtownodename;

  int m_cent_nclus_cut{0};

  // histos lists
  //  std::arrays have their indices backward, this is the old TH1 *cemc_hist_eta_phi[96][258];
  std::array<std::array<TH1 *, 258>, 96> cemc_hist_eta_phi{};
  std::array<TH1 *, 96> eta_hist{};
  TH2F *mass_eta{nullptr};
  TH3F *mass_eta_phi{nullptr};
  TH1F *h_totalClusters{nullptr};
  TH3F *pt1_ptpi0_alpha{nullptr};
  TH2F *fitp1_eta_phi2d{nullptr};
  TH1F *pairInvMassTotal{nullptr};

  TTree *_eventTree{nullptr};
  // TTree variables
  int _eventNumber{-1};
  int _nClusters{-1};
  float _clusterIDs[10000]{0};
  float _clusterEnergies[10000]{0};
  float _clusterPts[10000]{0};
  float _clusterEtas[10000]{0};
  float _clusterPhis[10000]{0};

  int maxTowerEta{-1};
  int maxTowerPhi{-1};

  int _maxTowerEtas[10000]{0};
  int _maxTowerPhis[10000]{0};

  float alphaCut{-1.};
  // TNtuple -> to store fit parameters

  /* TNtuple *nt_corrVals; */
  /* TF1 *fit_func; */
  /* TF1 *fit_result; */
  /* float fit_value_mean; */
  /* float corr_val; */

  TFile *f_temp{nullptr};

  int m_UseTowerInfo{0};  // 0 only old tower, 1 only new (TowerInfo based),
};

#endif  //   CALOEMCPI0TBT_CALOCALIBEMC_PI0_H
