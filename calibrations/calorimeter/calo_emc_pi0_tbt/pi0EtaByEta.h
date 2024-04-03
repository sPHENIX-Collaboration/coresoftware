#ifndef CALOANA_H__
#define CALOANA_H__

#include <fun4all/SubsysReco.h>

//#include <CLHEP/Vector/ThreeVector.h>  // for Hep3Vector
#include <array>
#include <string>  // for string
#include <vector>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TH2;
class TH1;
class TF1;
class TProfile2D;
class TH3;

namespace CLHEP
{
  class Hep3Vector;
}

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
  void set_timing_cut_width(const int t) { _range = t; }
  void set_vertex_cut(const float v) { _vz = v; }
  void apply_vertex_cut(bool Vtx_cut) { m_vtxCut = Vtx_cut; }

  TF1* fitHistogram(TH1* h);
  void fitEtaSlices(const std::string& infile, const std::string& outfile, const std::string& cdbFile);
  void set_use_pdc(bool state)
  {
    use_pdc = state;
    return;
  }

  void set_pt1BaseClusCut(float fac)
  {
    pt1BaseClusCut = fac;
    return;
  }
  void set_pt2BaseClusCut(float fac)
  {
    pt2BaseClusCut = fac;
    return;
  }
  void set_NclusDeptFac(float fac)
  {
    NclusDeptFac = fac;
    return;
  }
  void set_doMix(bool state)
  {
    doMix = state;
    return;
  }

  void set_massTargetHistFile(const std::string& file);

 protected:
  int Getpeaktime(TH1* h);
  std::string detector;
  std::string outfilename;

  float pt1BaseClusCut = 1.3;
  float pt2BaseClusCut = 0.7;
  float NclusDeptFac = 1.4;

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

  std::array<TH1*, 96> h_mass_eta_lt{};

  int _eventcounter{0};
  int _range{1};

  float _vz{0.};
  float target_pi0_mass{0.152};

  bool m_vtxCut{false};
  bool dynMaskClus{false};
  bool doMix{false};
  bool use_pdc{false};

  std::vector<std::vector<std::vector<CLHEP::Hep3Vector>>>* clusMix;
  TH1* h_nclus_bin{nullptr};
  const int NBinsClus = 10;
  TH1* h_vtx_bin{nullptr};
  int NBinsVtx = 30;
  TH1* h_event{nullptr};

  Fun4AllHistoManager* hm{nullptr};
  TFile* outfile{nullptr};
  TH2* h_emcal_mbd_correlation{nullptr};
  TH2* h_ohcal_mbd_correlation{nullptr};
  TH2* h_ihcal_mbd_correlation{nullptr};
  TH2* h_emcal_hcal_correlation{nullptr};
  TH2* h_emcal_zdc_correlation{nullptr};
  std::array<TH1*, 100> h_InvMass_Nclus{};

  TH1* h_InvMass{nullptr};
  TH1* h_InvMassMix{nullptr};

  TH1* h_target_mass{nullptr};
  bool use_h_target_mass{false};

  TH2* h_cemc_etaphi{nullptr};
  TH2* h_hcalin_etaphi{nullptr};
  TH2* h_hcalout_etaphi{nullptr};
  TH2* h_cemc_etaphi_wQA{nullptr};
  TH2* h_hcalin_etaphi_wQA{nullptr};
  TH2* h_hcalout_etaphi_wQA{nullptr};
  TH1* h_totalzdc_e{nullptr};
  TH3* h_pipT_Nclus_mass{nullptr};

  TProfile2D* h_cemc_etaphi_time{nullptr};
  TProfile2D* h_hcalin_etaphi_time{nullptr};
  TProfile2D* h_hcalout_etaphi_time{nullptr};

  TProfile2D* h_cemc_etaphi_badChi2{nullptr};
  TProfile2D* h_hcalin_etaphi_badChi2{nullptr};
  TProfile2D* h_hcalout_etaphi_badChi2{nullptr};

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

  TH1* h_clus_pt{nullptr};
  TH2* h_etaphi_clus{nullptr};

  TNtuple* g4hitntuple{nullptr};
  TNtuple* g4cellntuple{nullptr};
  TTree* towerntuple{nullptr};
  TNtuple* clusterntuple{nullptr};
  TH1* h_cemc_etaphi_noCalib{nullptr};

  TH1* h_pt1{nullptr};
  TH1* h_pt2{nullptr};
  TH1* h_nclusters{nullptr};
  TH1* h_emcal_e_eta{nullptr};
};

#endif
