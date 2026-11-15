#ifndef CALOVALID_CDBCALOQA_H
#define CALOVALID_CDBCALOQA_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <map>
#include <string>
#include <vector>

class CDBTTree;
class Fun4AllHistoManager;
class PHCompositeNode;
class TH1;
class TH2;
class TProfile;
class TProfile2D;
class TriggerAnalyzer;

// calorimeter, MBD, ZDC and sEPD histograms of one run for the CDB QA pages;
// in select mode it only drops the events not needed, before the reconstruction
class CDBCaloQA : public SubsysReco
{
 public:
  CDBCaloQA(const std::string &name = "CDBCaloQA", const std::string &outfile = "CDBCaloQA.root");
  ~CDBCaloQA() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_trigger_MB_bit(int bit) { m_mbbit_user = bit; }
  void set_MB_or_prescale(int ps) { m_mbor_prescale = ps; }
  void set_trigger_photon(const std::string &n) { m_trigname[kPhoton] = n; }
  void set_trigger_jet(const std::string &n) { m_trigname[kJet] = n; }

  void set_min_mb(int n) { m_minmb = n; }
  void set_min_zdc(int n) { m_minzdc = n; }
  void set_min_trig(int n) { m_mintrig = n; }
  void set_select_only(bool b) { m_select = b; }
  void set_centrality_vz_cut(float cut) { m_cent_vzcut = cut; }

  void set_tower_threshold(int layer, float e) { m_towthr.at(layer) = e; }
  void set_ns_per_sample(float ns) { m_ns_per_sample = ns; }
  void set_noZSXCalib_prefix(const std::string &p) { m_noxc_prefix = p; }
  void set_jet_node(const std::string &n) { m_jetnode = n; }
  void set_mc(bool b) { m_mc = b; }
  void set_species(const std::string &s) { m_species_user = s; }

  enum TrigCat
  {
    kMB = 0,
    kPhoton = 1,
    kJet = 2,
    kNTrig = 3
  };
  enum Layer
  {
    kCEMC = 0,
    kHCALIN = 1,
    kHCALOUT = 2,
    kNLayer = 3
  };

 private:
  void createHistos();
  void captureCdbUrls();
  void processMbd(PHCompositeNode *topNode, bool isMB);
  void processTowers(PHCompositeNode *topNode, int layer);
  void processClusters(PHCompositeNode *topNode, const std::array<bool, kNTrig> &fired, bool isMB);
  void processJets(PHCompositeNode *topNode, const std::array<bool, kNTrig> &fired);
  void processCentrality(PHCompositeNode *topNode, bool isMB);
  void processZdc(PHCompositeNode *topNode);
  void processSepd(PHCompositeNode *topNode, bool isMB);
  TH1 *book1(const std::string &name, const std::string &title, int nb, double lo, double hi);
  TH2 *book2(const std::string &name, const std::string &title, int nx, double xlo, double xhi, int ny, double ylo, double yhi);

  std::string m_outfile;
  std::string m_species{"pp"};
  std::string m_species_user;
  bool m_mc{false};
  int m_runnumber{0};
  bool m_booked{false};
  bool m_cdbcaptured{false};
  std::map<std::string, std::string> m_cdburls;

  static std::vector<int> zdcBits(TriggerAnalyzer *trigana);
  static std::vector<int> mbdBits(TriggerAnalyzer *trigana);
  static bool zdcEvent(TriggerAnalyzer *trigana, const std::vector<int> &zdcbits, const std::vector<int> &vetobits);

  bool m_select{false};
  int m_mbbit_user{-1};
  int m_minmb{0};
  int m_minzdc{0};
  int m_mintrig{0};
  int m_nevt{0};
  int m_nmb{0};
  int m_nzdc{0};
  std::array<int, kNTrig> m_ntrig{};
  bool m_zdcdecided{false};
  std::vector<int> m_zdcbits;
  std::vector<int> m_mbdbits;
  std::string m_zdcname;
  int m_mbbit{10};
  int m_mbbit_or{-1};
  int m_mbor_prescale{1000};
  std::string m_mbname;
  std::array<std::string, kNTrig> m_trigname{"", "Photon 4 GeV + MBD NS >= 1", "Jet 10 GeV + MBD NS >= 1"};
  std::array<float, kNLayer> m_towthr{0.1, 0.1, 0.1};
  float m_ns_per_sample{17.6};
  float m_noise_vzcut{100.};
  float m_cent_vzcut{1e9};
  bool m_noise_evt{false};
  std::string m_noxc_prefix{"TOWERINFO_NOZSXCALIB_"};
  std::string m_jetnode{"AntiKt_unsubtracted_r04"};

  float m_vtxz{0};
  float m_mbdt0{0};

  Fun4AllHistoManager *m_hm{nullptr};
  TriggerAnalyzer *m_trigana{nullptr};

  TH1 *h_evt{nullptr};
  TH1 *h_runnumber{nullptr};
  TH1 *h_vtxz{nullptr};

  TH2 *h2_mbd_q_pmt{nullptr};
  TH2 *h2_mbd_t_pmt{nullptr};
  TH1 *h_mbd_qavg_south{nullptr};
  TH1 *h_mbd_qavg_north{nullptr};
  TH1 *h_mbd_npmt_south{nullptr};
  TH1 *h_mbd_npmt_north{nullptr};
  TH1 *h_mbd_zvtx{nullptr};
  TH1 *h_mbd_t0{nullptr};

  std::array<TH1 *, kNLayer> h_tow_e{};
  std::array<TH1 *, kNLayer> h_tow_e_noZSXCalib{};
  std::array<TH1 *, kNLayer> h_tow_esum{};
  std::array<TH1 *, kNLayer> h_tow_occ{};
  std::array<TProfile *, kNLayer> p_tow_occ_mean{};
  std::array<TProfile *, kNLayer> p_tsc_e{};
  std::array<TProfile *, kNLayer> p_tsc_occ{};
  std::array<float, kNLayer> m_tsclo{0, 0.3, 1.0};
  std::array<float, kNLayer> m_tschi{0, 1.0, 2.5};
  std::array<TH1 *, kNLayer> h_tow_nonzs_frac{};
  std::array<TH1 *, kNLayer> h_tow_time{};
  std::array<TH1 *, kNLayer> h_tow_time_mbdt0{};
  TH1 *h_noise_nevt{nullptr};
  std::array<TH1 *, kNLayer> h_noise_e_fine{};
  std::array<TH2 *, kNLayer> h2_noise_e_fine_ieta{};

  TH1 *h_cent{nullptr};
  TH1 *h_cent_mbd_qsum{nullptr};

  CDBTTree *m_zdccalib{nullptr};
  TH1 *h_zdc_e_south{nullptr};
  TH1 *h_zdc_e_north{nullptr};
  std::array<TH1 *, 6> h_zdc_mod_e{};
  TH2 *h2_zdc_e_ns{nullptr};
  TH2 *h2_zdc_e_ns_mb{nullptr};
  TH1 *h_zdc_e_south_mb{nullptr};
  TH1 *h_zdc_e_north_mb{nullptr};
  std::array<TH2 *, 2> h2_zdc_e_ohcal{};
  std::array<TProfile2D *, 2> p2_zdc_t0_ohcal{};
  std::array<TH2 *, 2> h2_zdc_e_t0{};

  float m_sepd_hitthr{0.5};
  TH2 *h2_sepd_q_ch{nullptr};
  TH1 *h_sepd_qsum_south{nullptr};
  TH1 *h_sepd_qsum_north{nullptr};
  TH2 *h2_sepd_qsum_ns{nullptr};
  TH1 *h_sepd_nhit_south{nullptr};
  TH1 *h_sepd_nhit_north{nullptr};

  std::array<TH2 *, kNTrig> h2_clus_e_eta{};
  TH1 *h_pi0_mass{nullptr};
  std::array<TH2 *, kNTrig> h2_jet_pt_eta{};
};

#endif
