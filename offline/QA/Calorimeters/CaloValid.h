#ifndef CALOVALID_CALOVALID_H
#define CALOVALID_CALOVALID_H

#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TH1;
class TH2;
class TProfile2D;

class CaloValid : public SubsysReco
{
 public:
  //! constructor
  CaloValid(const std::string& name = "CaloValid");// const std::string &filename = "testQA.root"); //int nevents = 100);

  //! destructor
  virtual ~CaloValid();

  //! full initialization
  int Init(PHCompositeNode*);

  //! event processing method
  int process_event(PHCompositeNode*);

  //! end of run method
  int End(PHCompositeNode*);

  int process_g4hits(PHCompositeNode*);
  int process_g4cells(PHCompositeNode*);
  int process_towers(PHCompositeNode*);
  int process_clusters(PHCompositeNode*);

  void Detector(const std::string& name) { detector = name; }
  void set_timing_cut_width(const int& t) { _range = t; }

  void set_debug(bool debug) { m_debug = debug; }
  TH2* LogYHist2D(const std::string& name, const std::string& title, int, double, double, int, double, double);

 private:
  int Getpeaktime(TH1* h);
  void createHistos();
  void MirrorHistogram(TH1* histogram);
  std::string getHistoPrefix() const;
  bool m_debug{0};
  std::string detector;
  TFile* OutputNtupleFile;
  std::string m_outputFileName;
  std::string OutputFileName;

  TH1* h_cemc_channel_pedestal[128*192];
  TH1* h_ihcal_channel_pedestal[32*48];
  TH1* h_ohcal_channel_pedestal[32*48];

  TH1* h_cemc_channel_energy[128*192];
  TH1* h_ihcal_channel_energy[32*48];
  TH1* h_ohcal_channel_energy[32*48];

  //TProfile2D* h_cemc_etaphi_pedRMS{nullptr};
 
  //TProfile2D* h_cemc_etaphi_pedRMS{nullptr};
  /*
  TH2* h_emcal_mbd_correlation{nullptr};
  TH2* h_ohcal_mbd_correlation{nullptr};
  TH2* h_ihcal_mbd_correlation{nullptr};
  TH2* h_emcal_hcal_correlation{nullptr};
  TH2* h_emcal_zdc_correlation{nullptr};

  TH1* h_InvMass{nullptr};

  TH2* h_cemc_etaphi{nullptr};
  TH2* h_hcalin_etaphi{nullptr};
  TH2* h_hcalout_etaphi{nullptr};
  TH2* h_cemc_etaphi_wQA{nullptr};
  TH2* h_hcalin_etaphi_wQA{nullptr};
  TH2* h_hcalout_etaphi_wQA{nullptr};
  TH1* h_totalzdc_e{nullptr};

  TProfile2D* h_cemc_etaphi_time{nullptr};
  TProfile2D* h_hcalin_etaphi_time{nullptr};
  TProfile2D* h_hcalout_etaphi_time{nullptr};

  TH2* h_cemc_e_chi2{nullptr};
  TH2* h_ohcal_e_chi2{nullptr};
  TH2* h_ihcal_e_chi2{nullptr};

  TProfile2D* h_cemc_etaphi_badChi2{nullptr};
  TProfile2D* h_hcalin_etaphi_badChi2{nullptr};
  TProfile2D* h_hcalout_etaphi_badChi2{nullptr};

  TProfile2D* h_cemc_etaphi_fracHitADC{nullptr};
  TProfile2D* h_hcalin_etaphi_fracHitADC{nullptr};
  TProfile2D* h_hcalout_etaphi_fracHitADC{nullptr};

  TH1* hzdctime{nullptr};
  TH1* hmbdtime{nullptr};
  TH1* hemcaltime{nullptr};
  TH1* hihcaltime{nullptr};
  TH1* hohcaltime{nullptr};

  TH1* hzdctime_cut{nullptr};
  TH1* hmbdtime_cut{nullptr};
  TH1* hemcaltime_cut{nullptr};
  TH1* hihcaltime_cut{nullptr};
  TH1* hohcaltime_cut{nullptr};

  TH1* hvtx_z_raw{nullptr};
  TH1* hvtx_z_cut{nullptr};

  TH1* hzdcSouthraw{nullptr};
  TH1* hzdcNorthraw{nullptr};
  TH1* hzdcSouthcalib{nullptr};
  TH1* hzdcNorthcalib{nullptr};
  TH1* h_ihcal_status{nullptr};
  TH1* h_ohcal_status{nullptr};
  TH1* h_cemc_status{nullptr};

  TH1* h_clusE{nullptr};
  TH2* h_etaphi_clus{nullptr};

  TTree* towerntuple{nullptr};
  TNtuple* clusterntuple{nullptr};
*/

  int _eventcounter{0};
  int _range{1};
};

#endif
