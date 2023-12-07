#ifndef CALOANA_H__
#define CALOANA_H__

#include <fun4all/SubsysReco.h>
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
class TLorentzVector;

class pi0EtaByEta : public SubsysReco
{
 public:
  //! constructor
  pi0EtaByEta(const std::string& name = "pi0EtaByEta", const std::string& fname = "MyNtuple.root");

  //! destructor
  virtual ~pi0EtaByEta();

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

  std::pair<double, double> fitHistogram(TH1F* h);
  void fitEtaSlices(std::string infile, std::string outfile, std::string cdbFile);

 protected:
  std::string detector;
  std::string outfilename;
  int Getpeaktime(TH1* h);
  Fun4AllHistoManager* hm = nullptr;
  TFile* outfile = nullptr;
  TH2F* h_emcal_mbd_correlation = nullptr;
  TH2F* h_ohcal_mbd_correlation = nullptr;
  TH2F* h_ihcal_mbd_correlation = nullptr;
  TH2F* h_emcal_hcal_correlation = nullptr;
  TH2F* h_emcal_zdc_correlation = nullptr;

  TH1F* h_InvMass = nullptr;
  TH1F* h_InvMassMix = nullptr;

  TH2F* h_cemc_etaphi = nullptr;
  TH2F* h_hcalin_etaphi = nullptr;
  TH2F* h_hcalout_etaphi = nullptr;
  TH2F* h_cemc_etaphi_wQA = nullptr;
  TH2F* h_hcalin_etaphi_wQA = nullptr;
  TH2F* h_hcalout_etaphi_wQA = nullptr;
  TH1* h_totalzdc_e;

  TProfile2D* h_cemc_etaphi_time = nullptr;
  TProfile2D* h_hcalin_etaphi_time = nullptr;
  TProfile2D* h_hcalout_etaphi_time = nullptr;

  TProfile2D* h_cemc_etaphi_badChi2 = nullptr;
  TProfile2D* h_hcalin_etaphi_badChi2 = nullptr;
  TProfile2D* h_hcalout_etaphi_badChi2 = nullptr;

  TH1* hzdctime;
  TH1* hmbdtime;
  TH1* hemcaltime;
  TH1* hihcaltime;
  TH1* hohcaltime;

  TH1* hzdctime_cut;
  TH1* hmbdtime_cut;
  TH1* hemcaltime_cut;
  TH1* hihcaltime_cut;
  TH1* hohcaltime_cut;

  TH1* hvtx_z_raw;
  TH1* hvtx_z_cut;

  TH1* hzdcSouthraw;
  TH1* hzdcNorthraw;
  TH1* hzdcSouthcalib;
  TH1* hzdcNorthcalib;

  TH1F* h_clusE;
  TH2F* h_etaphi_clus;

  TNtuple* g4hitntuple = nullptr;
  TNtuple* g4cellntuple = nullptr;
  TTree* towerntuple = nullptr;
  TNtuple* clusterntuple = nullptr;
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
  int _eventcounter;
  int _range = 1;
  float _vz = 0.;
  bool m_vtxCut = false;
  bool dynMaskClus = true;

  TH1F* h_pt1;
  TH1F* h_pt2;
  TH1F* h_nclusters;
  TH1F* h_mass_eta_lt[96];
  TH1F* h_emcal_e_eta;

  float target_pi0_mass = 0.145;
};

#endif
