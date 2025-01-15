#ifndef CALOVALID_CALOVALID_H
#define CALOVALID_CALOVALID_H

#include <calotrigger/TriggerAnalyzer.h>

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

// Forward declarations
class PHCompositeNode;
class TH1;
class TH2;
class TH3;
class TProfile2D;
class TProfile;

class CaloValid : public SubsysReco
{
 public:
  //! constructor
  CaloValid(const std::string& name = "CaloValid");  // const std::string &filename = "testQA.root"); //int nevents = 100);

  //! destructor
  ~CaloValid() override;

  //! full initialization
  int Init(PHCompositeNode*) override;

  //! event processing method
  int process_event(PHCompositeNode*) override;

  //! end of run method
  int End(PHCompositeNode*) override;

  int process_g4hits(PHCompositeNode*);
  int process_g4cells(PHCompositeNode*);
  int process_towers(PHCompositeNode*);
  int process_clusters(PHCompositeNode*);

  void set_timing_cut_width(const int& t) { _range = t; }

  void set_debug(bool debug) { m_debug = debug; }
  TH2* LogYHist2D(const std::string& name, const std::string& title, int, double, double, int, double, double);

 private:
  void createHistos();
  void MirrorHistogram(TH1* histogram);
  std::string getHistoPrefix() const;

  TriggerAnalyzer* trigAna{nullptr};
  TH3* h_pi0_trigIB_mass {nullptr};
  std::vector<int> triggerIndices {10, 28, 29, 30, 31}; //MBD NS>=1, Photon Triggers
    
  TH1* h_cemc_channel_pedestal[128 * 192]{nullptr};
  TH1* h_ihcal_channel_pedestal[32 * 48]{nullptr};
  TH1* h_ohcal_channel_pedestal[32 * 48]{nullptr};

  TH1* h_cemc_channel_energy[128 * 192]{nullptr};
  TH1* h_ihcal_channel_energy[32 * 48]{nullptr};
  TH1* h_ohcal_channel_energy[32 * 48]{nullptr};
  TH2* h_emcal_mbd_correlation{nullptr};
  TH1* h_mbd_hits{nullptr};
  TH2* h_ohcal_mbd_correlation{nullptr};
  TH2* h_ihcal_mbd_correlation{nullptr};
  TH2* h_emcal_hcal_correlation{nullptr};
  TH2* h_cemc_etaphi{nullptr};
  TH2* h_ihcal_etaphi{nullptr};
  TH2* h_ohcal_etaphi{nullptr};
  TH2* h_cemc_etaphi_wQA{nullptr};
  TH2* h_ihcal_etaphi_wQA{nullptr};
  TH2* h_ohcal_etaphi_wQA{nullptr};
  TH1* h_ihcal_status{nullptr};
  TH1* h_ohcal_status{nullptr};
  TH1* h_cemc_status{nullptr};
  TH2* h_cemc_e_chi2{nullptr};
  TH2* h_ihcal_e_chi2{nullptr};
  TH2* h_ohcal_e_chi2{nullptr};
  TProfile2D* h_cemc_etaphi_time{nullptr};
  TProfile2D* h_ihcal_etaphi_time{nullptr};
  TProfile2D* h_ohcal_etaphi_time{nullptr};
  TProfile2D* h_cemc_etaphi_fracHitADC{nullptr};
  TProfile2D* h_cemc_etaphi_fracHit{nullptr};
  TProfile2D* h_ihcal_etaphi_fracHitADC{nullptr};
  TProfile2D* h_ohcal_etaphi_fracHitADC{nullptr};
  TProfile2D* h_cemc_etaphi_pedRMS{nullptr};
  TProfile2D* h_ihcal_etaphi_pedRMS{nullptr};
  TProfile2D* h_ohcal_etaphi_pedRMS{nullptr};
  TProfile2D* h_cemc_etaphi_ZSpedRMS{nullptr};
  TProfile2D* h_ohcal_etaphi_ZSpedRMS{nullptr};
  TProfile2D* h_ihcal_etaphi_ZSpedRMS{nullptr};
  TProfile2D* h_cemc_etaphi_badChi2{nullptr};
  TProfile2D* h_ohcal_etaphi_badChi2{nullptr};
  TProfile2D* h_ihcal_etaphi_badChi2{nullptr};
  TH1* h_InvMass{nullptr};
  TH1* h_channel_pedestal_0{nullptr};
  TH1* h_vtx_z_raw{nullptr};
  TH1* h_vtx_z_cut{nullptr};
  TH1* h_emcaltime_cut{nullptr};
  TH1* h_ihcaltime_cut{nullptr};
  TH1* h_ohcaltime_cut{nullptr};
  TH1* h_emcaltime{nullptr};
  TH1* h_ihcaltime{nullptr};
  TH1* h_ohcaltime{nullptr};
  TH1* h_emcal_tower_e{nullptr};
  TH2* h_etaphi_clus{nullptr};
  TH1* h_clusE{nullptr};

  TProfile2D* h_cemc_etaphi_time_raw{nullptr};
  TProfile2D* h_ohcal_etaphi_time_raw{nullptr};
  TProfile2D* h_ihcal_etaphi_time_raw{nullptr};

  // Trigger histos
  TH1 *h_triggerVec{nullptr};
  TH2 *h_edist[64] = {nullptr};
  TH1 *h_ldClus_trig[64] = {nullptr};
  TProfile *pr_evtNum_ldClus_trig[64] = {nullptr};
  TProfile *pr_rejection[64] = {nullptr};
  TProfile *pr_livetime[64] = {nullptr};
  TProfile *pr_ldClus_trig{nullptr};
  std::vector<int> trigOfInterest = {3,10,11,21,22,23,25,26,27};

  int _eventcounter{0};
  int _range{1};

  bool m_debug{false};

  std::string m_outputFileName;
  std::string OutputFileName;
};

#endif

