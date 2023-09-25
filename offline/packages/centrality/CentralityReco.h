#ifndef CENTRALITY_CENTRALITYRECO_H
#define CENTRALITY_CENTRALITYRECO_H

#include <fun4all/SubsysReco.h>

#include <limits>

// Forward declarations
class CentralityInfo;
class Fun4AllHistoManager;
class PHCompositeNode;
class TProfile;
class TFile;
class TNtuple;
class TTree;
class TowerInfo;
class TowerInfoContainer;
class TH1;
class TH2;
class BbcVertex;
class BbcVertexMapv1;

class CentralityReco : public SubsysReco
{
 public:
  //! constructor
  explicit CentralityReco(const std::string &name = "CentralityReco", const std::string &hist_name = "QA_CentralityReco.root", const std::string &tree_name = "centrality.root");

  //! destructor
  virtual ~CentralityReco();

  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  void CreateNodes(PHCompositeNode *);
  void PrintCentiles();
  int CheckZDC();
  int GetMBDVertexAndCharge();
  int FillCentralityInfo();
  void FillHistograms();
  int FillVars();
  int GetNodes(PHCompositeNode *);
  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  enum OperationMode
  {
    QA = 0,
    Performance = 1
  };

  void SetOperationMode(int imode) { _op_mode = imode; }
  void ResetVars();

 private:
  Fun4AllHistoManager *hm = nullptr;

  TFile *outfile = nullptr;
  TTree *ttree = nullptr;
  TowerInfo *_tmp_tower = nullptr;
  BbcVertexMapv1 *_bbc_vertex_map = nullptr;
  CentralityInfo *_central = nullptr;
  TowerInfoContainer *_towers_mbd = nullptr;
  TowerInfoContainer *_towers_zdc = nullptr;

  // Histograms
  //  mbd
  TH1 *h_mbd_vertex = nullptr;
  TH1 *h_mbd_vertex_w_zdc_cut = nullptr;

  TH1 *h_mbd_charge_ns = nullptr;
  TH1 *h_mbd_charge_n = nullptr;
  TH1 *h_mbd_charge_s = nullptr;
  TH1 *h_mbd_ring_charge_sum_n[3]{};
  TH1 *h_mbd_ring_charge_sum_s[3]{};

  TH1 *h_mbd_charge_ns_w_zdc_cut = nullptr;
  TH1 *h_mbd_charge_n_w_zdc_cut = nullptr;
  TH1 *h_mbd_charge_s_w_zdc_cut = nullptr;
  TH1 *h_mbd_ring_charge_sum_n_w_zdc_cut[3]{};
  TH1 *h_mbd_ring_charge_sum_s_w_zdc_cut[3]{};

  TH1 *h_mbd_charge_ns_w_zdc_cut_w_mbd_cut = nullptr;
  TH1 *h_mbd_charge_n_w_zdc_cut_w_mbd_cut = nullptr;
  TH1 *h_mbd_charge_s_w_zdc_cut_w_mbd_cut = nullptr;
  TH1 *h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[3]{};
  TH1 *h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[3]{};

  TH1 *h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[5]{};
  TH1 *h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[5]{};
  TH1 *h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[5]{};
  TH1 *h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[3][5]{};
  TH1 *h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[3][5]{};

  // zdc
  TH1 *h_zdc_energy_ns = nullptr;
  TH1 *h_zdc_energy_n = nullptr;
  TH1 *h_zdc_energy_s = nullptr;
  /* TH1 *h_zdc_energy_single_n[3] {}; */
  /* TH1 *h_zdc_energy_single_s[3] {}; */

  // correlations
  TH2 *h_zdc_mbd_corr_ns = nullptr;
  TH2 *h_zdc_mbd_corr_n = nullptr;
  TH2 *h_zdc_mbd_corr_s = nullptr;

  TH2 *h_zdc_mbd_corr_ns_w_zdc_cut = nullptr;
  TH2 *h_zdc_mbd_corr_n_w_zdc_cut = nullptr;
  TH2 *h_zdc_mbd_corr_s_w_zdc_cut = nullptr;

  TH2 *h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut = nullptr;
  TH2 *h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut = nullptr;
  TH2 *h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut = nullptr;

  TH2 *h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[5]{};
  TH2 *h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[5]{};
  TH2 *h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[5]{};

  std::string _hist_filename;
  std::string _tree_filename;

  bool _zdc_check = false;

  unsigned int _key = std::numeric_limits<unsigned int>::max();

  int _op_mode = 0;
  int _isMinBias = std::numeric_limits<int>::max();
  int _quality = std::numeric_limits<int>::max();
  int _side = std::numeric_limits<int>::max();
  int _channel = std::numeric_limits<int>::max();
  int _type = std::numeric_limits<int>::max();

  const int mbd_ring_index[64] =
      {2, 2, 2, 1, 1, 2, 1, 0,
       0, 2, 1, 0, 2, 1, 0, 1,
       0, 2, 1, 0, 2, 1, 0, 2,
       1, 0, 0, 2, 1, 1, 2, 2,
       2, 2, 2, 1, 1, 2, 1, 0,
       0, 2, 1, 0, 2, 1, 0, 1,
       0, 2, 1, 0, 2, 1, 0, 2,
       1, 0, 0, 2, 1, 1, 2, 2};

  int _tubes_hit[2]{};
  int _tdc[2]{};
  int m_mbd_side[128]{};
  int m_mbd_channel[128]{};

  float _mbd_charge_threshold = 0.3;
  float _mbd_time_threshold = 15.0;
  float _zdc_energy_threshold = 1.;
  float _z_vertex = std::numeric_limits<float>::signaling_NaN();
  float _offset = std::numeric_limits<float>::signaling_NaN();
  float _energy = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum_n = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum_s = std::numeric_limits<float>::signaling_NaN();
  float _zdc_energy_sum = std::numeric_limits<float>::signaling_NaN();
  float _zdc_energy_sum_n = std::numeric_limits<float>::signaling_NaN();
  float _zdc_energy_sum_s = std::numeric_limits<float>::signaling_NaN();

  float _mbd_vertex_cuts[5]{};
  float _centrality_map[20]{};
  float _zdc_gain_factors[6]{};
  float gaincorr[128]{};
  float tq_t0_offsets[128]{};
  float _mbd_ring_charge_sum_n[3]{};
  float _mbd_ring_charge_sum_s[3]{};

  float m_mbd_charge[128]{};
  float m_mbd_time[128]{};
  float m_mbd_charge_raw[128]{};
  float m_mbd_time_raw[128]{};

  float m_zdc_energy_low[6]{};
  float m_zdc_energy_high[6]{};
  float m_zdc_sum_low[2]{};
  float m_zdc_sum_high[2]{};
};

#endif
