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
class TH2F;
class TH1F;
class TH1;
class TProfile2D;

class CaloValid : public SubsysReco
{
 public:
  //! constructor
  CaloValid(const std::string& name = "CaloValid", const std::string& fname = "MyNtuple.root");

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
  void set_vertex_cut(const float& v) { _vz = v; }
  void apply_vertex_cut(bool Vtx_cut) { m_vtxCut = Vtx_cut; }

  void set_debug(bool debug) { m_debug = debug; }
  TH2F* LogYHist2D(const std::string& name, const std::string& title, int, double, double, int, double, double);

 private:
  int Getpeaktime(TH1* h);

  bool m_debug{0};
  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager* hm{nullptr};
  TFile* outfile{nullptr};
  TH2F* h_emcal_mbd_correlation{nullptr};
  TH2F* h_ohcal_mbd_correlation{nullptr};
  TH2F* h_ihcal_mbd_correlation{nullptr};
  TH2F* h_emcal_hcal_correlation{nullptr};
  TH2F* h_emcal_zdc_correlation{nullptr};

  TH1F* h_InvMass{nullptr};

  TH2F* h_cemc_etaphi{nullptr};
  TH2F* h_hcalin_etaphi{nullptr};
  TH2F* h_hcalout_etaphi{nullptr};
  TH2F* h_cemc_etaphi_wQA{nullptr};
  TH2F* h_hcalin_etaphi_wQA{nullptr};
  TH2F* h_hcalout_etaphi_wQA{nullptr};
  TH1* h_totalzdc_e{nullptr};

  TProfile2D* h_cemc_etaphi_time{nullptr};
  TProfile2D* h_hcalin_etaphi_time{nullptr};
  TProfile2D* h_hcalout_etaphi_time{nullptr};

  TH2F* h_cemc_e_chi2{nullptr};
  TH2F* h_ohcal_e_chi2{nullptr};
  TH2F* h_ihcal_e_chi2{nullptr};

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

  TH1F* h_clusE{nullptr};
  TH2F* h_etaphi_clus{nullptr};

  TTree* towerntuple{nullptr};
  TNtuple* clusterntuple{nullptr};

  std::vector<float> m_energy;
  std::vector<int> m_etabin;
  std::vector<int> m_phibin;
  std::vector<int> m_time;

  std::vector<float> m_hcalin_energy;
  std::vector<int> m_hcalin_etabin;
  std::vector<int> m_hcalin_phibin;
  std::vector<int> m_hcalin_time;

  std::vector<float> m_hcalout_energy;
  std::vector<int> m_hcalout_etabin;
  std::vector<int> m_hcalout_phibin;
  std::vector<int> m_hcalout_time;

  std::vector<float> m_zdc_energy;
  std::vector<int> m_zdc_index;
  std::vector<int> m_zdc_side;

  std::vector<float> m_bbc_energy;
  std::vector<int> m_bbc_type;
  std::vector<int> m_bbc_side;
  int _eventcounter{0};
  int _range{1};
  float _vz{0.};
  bool m_vtxCut{false};
  bool dynMaskClus{true};
};

#endif
