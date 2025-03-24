// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOEMCNOISYTOWER_EMCNOISYTOWERFINDER_H
#define CALOEMCNOISYTOWER_EMCNOISYTOWERFINDER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class TFile;
class TH2;
class TProfile2D;
class CDBTTree;

class emcNoisyTowerFinder : public SubsysReco
{
 public:
  explicit emcNoisyTowerFinder(const std::string& name = "emcNoisyTowerFinder", const std::string& outputName = "emcNoisyTowerFinder.root");

  ~emcNoisyTowerFinder() override = default;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;
  void FillHistograms(const int runnumber, const int segment);
  void CalculateCutOffs(const int runnumber);
  void WriteCDBTree(const int runnumber);
  void FindHot(std::string& infilename, std::string& outfilename, const std::string& inHist = "h_hits_eta_phi_gev");

  void set_energy_threshold_adc(float val) { energy_threshold_adc = val; }
  void set_energy_threshold_gev(float val) { energy_threshold_gev = val; }
  void set_sigma_bad_thresh(float val) { sigma_bad_thresh = val; }
  void set_ihcal()
  {
    Nphi = 64;
    Neta = 24;
    m_caloName = "HCALIN";
    return;
  }
  void set_ohcal()
  {
    Nphi = 64;
    Neta = 24;
    m_caloName = "HCALOUT";
    return;
  }
  static float findMedian(const std::vector<float>& arr);

 private:
  TFile* out{nullptr};
  TFile* foutput{nullptr};

  CDBTTree* cdbttree{nullptr};

  TH2* h_hits_eta_phi_adc{nullptr};
  TProfile2D* pr_hits_eta_phi_adc{nullptr};
  TH2* h_hits_eta_phi_gev{nullptr};
  TProfile2D* pr_hits_eta_phi_gev{nullptr};

  int Neta{96};
  int Nphi{256};

  float energy_threshold_adc{200};
  float energy_threshold_gev{0.25};
  float percent_cold_thresh{0.1};
  float sigma_bad_thresh{5};

  std::string Outfile{"commissioning.root"};

  std::string m_fieldname{"Femc_datadriven_qm1_correction"};
  std::string m_caloName{"CEMC"};
};

#endif
