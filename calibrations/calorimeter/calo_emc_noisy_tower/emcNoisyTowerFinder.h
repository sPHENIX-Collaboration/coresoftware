// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOEMCNOISYTOWER_EMCNOISYTOWERFINDER_H
#define CALOEMCNOISYTOWER_EMCNOISYTOWERFINDER_H

#include <fun4all/SubsysReco.h>
// #include <cdbobjects/CDBTTree.h>

#include <array>
#include <string>
#include <vector>

class Fun4AllHistoManager;
class PHCompositeNode;
class RawCluster;
class TFile;
class TH2F;
class TowerInfoContainer;
class TProfile2D;
class CDBTTree;

class emcNoisyTowerFinder : public SubsysReco
{
 public:
  explicit emcNoisyTowerFinder(const std::string& name = "emcNoisyTowerFinder", const std::string& outputName = "emcNoisyTowerFinder.root");

  ~emcNoisyTowerFinder() override;
  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int ResetEvent(PHCompositeNode* topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode* topNode) override;
  int Reset(PHCompositeNode* /*topNode*/) override;
  void Print(const std::string& what = "ALL") const override;
  void FillHistograms(const int runnumber, const int segment);
  void CalculateCutOffs(const int runnumber);
  void WriteCDBTree(const int runnumber);
  void FindHot(std::string& infilename, std::string& outfilename, const std::string& inHist = "h_hits_eta_phi_gev");

  void set_energy_threshold_adc(float val) { energy_threshold_adc = val; }
  void set_energy_threshold_gev(float val) { energy_threshold_gev = val; }
  void set_sigma_bad_thresh(float val) { sigma_bad_thresh = val; }
  void set_hcal(){Nphi = 64; Neta = 24; return;}

 private:
  TFile* out{nullptr};

  CDBTTree* cdbttree{nullptr};

  std::string Outfile{"commissioning.root"};

  std::string m_fieldname = "Femc_datadriven_qm1_correction";

  TH2F* h_hits_eta_phi_adc{nullptr};
  TProfile2D* pr_hits_eta_phi_adc{nullptr};
  TH2F* h_hits_eta_phi_gev{nullptr};
  TProfile2D* pr_hits_eta_phi_gev{nullptr};

  int Nphi = 256;
  int Neta = 96;

  float energy_threshold_adc = 200;
  float energy_threshold_gev = 0.25;

  TFile* foutput{nullptr};

  float sigma_bad_thresh = 5;
  float percent_cold_thresh = 0.1;
};

#endif
