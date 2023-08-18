#ifndef CENTRALITYRECO_H__
#define CENTRALITYRECO_H__

#include "CentralityInfov1.h"
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <calobase/TowerInfoContainerv1.h>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TProfile;
class TFile;
class TNtuple;
class TTree;
class TowerInfoContainerv1;
class TH1;
class TH2;

class CentralityReco : public SubsysReco
{
 public:
  //! constructor
  CentralityReco(const std::string &name = "CentralityReco", const std::string &hist_name = "QA_CentralityReco.root", const std::string &tree_name = "centrality.root");

  //! destructor
  virtual ~CentralityReco();

  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  void CreateNodes(PHCompositeNode*);
  void PrintCentiles();
  int CheckZDC();
  int GetMBDVertexAndCharge();
  void FillHistograms();
  int FillVars();
  int GetNodes(PHCompositeNode *);
  void SetVerbosity(int verbose) {_verbose = verbose;}
  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);
  
  enum OperationMode
  {
    QA = 0,
    Performance = 1
  };

  void SetOperationMode(int imode) { _op_mode = imode;}
  void ResetVars();
 protected:

  Fun4AllHistoManager *hm = nullptr;
  std::string _hist_filename;

  TFile *outfile;
  TTree *ttree;
  std::string _tree_filename;

  int mbd_ring_index[64] = 
    {2,2,2,1,1,2,1,0,
     0,2,1,0,2,1,0,1,
     0,2,1,0,2,1,0,2,
     1,0,0,2,1,1,2,2,
     2,2,2,1,1,2,1,0,
     0,2,1,0,2,1,0,1,
     0,2,1,0,2,1,0,2,
     1,0,0,2,1,1,2,2};

  int _verbose;
  int _op_mode;
  float _mbd_charge_threshold;
  float _zdc_energy_threshold;
  float _mbd_vertex_cuts[5] {};
  
  std::vector<float> _centiles;

  float _z_vertex;
  int _isMinBias;
  int _quality;
  bool _zdc_check;
  int _tubes_hit[2] {};
  float _offset;

  TowerInfo *_tmp_tower;
  float _energy;
  unsigned int _key;
  int _side, _channel, _type;
  float _zdc_gain_factors[6];
  
  float gaincorr[128] {};
  float tq_t0_offsets[128] {};

  CentralityInfov1 *_central;
  TowerInfoContainerv1 *_towers_mbd;
  TowerInfoContainerv1 *_towers_zdc;
  
  float _mbd_charge_sum;
  float _mbd_charge_sum_n;
  float _mbd_charge_sum_s;
  float _mbd_ring_charge_sum_n[3];
  float _mbd_ring_charge_sum_s[3];

  float _zdc_energy_sum;
  float _zdc_energy_sum_n;
  float _zdc_energy_sum_s;

  float m_mbd_charge[128];
  float m_mbd_time[128];
  float m_mbd_charge_raw[128];
  float m_mbd_time_raw[128];

  int m_mbd_side[128];
  int m_mbd_channel[128];

  float m_zdc_energy_low[6];
  float m_zdc_energy_high[6];
  float m_zdc_sum_low[2];
  float m_zdc_sum_high[2];


  //Histograms
  // mbd
  TH1D *h_mbd_vertex;
  TH1D *h_mbd_vertex_w_zdc_cut;

  TH1D *h_mbd_charge_ns;
  TH1D *h_mbd_charge_n;
  TH1D *h_mbd_charge_s;
  TH1D *h_mbd_ring_charge_sum_n[3];
  TH1D *h_mbd_ring_charge_sum_s[3];

  TH1D *h_mbd_charge_ns_w_zdc_cut;
  TH1D *h_mbd_charge_n_w_zdc_cut;
  TH1D *h_mbd_charge_s_w_zdc_cut;
  TH1D *h_mbd_ring_charge_sum_n_w_zdc_cut[3];
  TH1D *h_mbd_ring_charge_sum_s_w_zdc_cut[3];

  TH1D *h_mbd_charge_ns_w_zdc_cut_w_mbd_cut;
  TH1D *h_mbd_charge_n_w_zdc_cut_w_mbd_cut;
  TH1D *h_mbd_charge_s_w_zdc_cut_w_mbd_cut;
  TH1D *h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut[3];
  TH1D *h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut[3];

  TH1D *h_mbd_charge_ns_w_zdc_cut_w_mbd_cut_and_vertex[5];
  TH1D *h_mbd_charge_n_w_zdc_cut_w_mbd_cut_and_vertex[5];
  TH1D *h_mbd_charge_s_w_zdc_cut_w_mbd_cut_and_vertex[5];
  TH1D *h_mbd_ring_charge_sum_n_w_zdc_cut_w_mbd_cut_and_vertex[3][5];
  TH1D *h_mbd_ring_charge_sum_s_w_zdc_cut_w_mbd_cut_and_vertex[3][5];

  // zdc
  TH1D *h_zdc_energy_ns;
  TH1D *h_zdc_energy_n;
  TH1D *h_zdc_energy_s;
  TH1D *h_zdc_energy_single_n[3];
  TH1D *h_zdc_energy_single_s[3];

  // correlations
  TH2D *h_zdc_mbd_corr_ns;
  TH2D *h_zdc_mbd_corr_n;
  TH2D *h_zdc_mbd_corr_s;

  TH2D *h_zdc_mbd_corr_ns_w_zdc_cut;
  TH2D *h_zdc_mbd_corr_n_w_zdc_cut;
  TH2D *h_zdc_mbd_corr_s_w_zdc_cut;

  TH2D *h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut;
  TH2D *h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut;
  TH2D *h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut;

  TH2D *h_zdc_mbd_corr_ns_w_zdc_cut_w_mbd_cut_and_vertex[5];
  TH2D *h_zdc_mbd_corr_n_w_zdc_cut_w_mbd_cut_and_vertex[5];
  TH2D *h_zdc_mbd_corr_s_w_zdc_cut_w_mbd_cut_and_vertex[5];

};

#endif
