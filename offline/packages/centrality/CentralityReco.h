#ifndef CENTRALITY_CENTRALITYRECO_H
#define CENTRALITY_CENTRALITYRECO_H

#include <fun4all/SubsysReco.h>

#include <limits>

// Forward declarations
class CentralityInfo;
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class TTree;
class TowerInfo;
class TowerInfoContainer;
class BbcPmtInfoV1;
class BbcPmtInfoContainerV1;
class BbcVertex;
class BbcVertexMapv1;

class CentralityReco : public SubsysReco
{
 public:
  //! constructor
  explicit CentralityReco(const std::string &name = "CentralityReco");

  //! destructor
  virtual ~CentralityReco();

  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  void CreateNodes(PHCompositeNode *);
  void PrintCentiles();
  int CheckZDC();
  int GetMBDCharge();
  int FillCentralityInfo();
  int FillVars();
  int GetNodes(PHCompositeNode *);
  //! event processing method
  int process_event(PHCompositeNode *);
  
  //! end of run method
  int End(PHCompositeNode *);

  void ResetVars();

  void useZDC(int use_zdc) { _use_ZDC = use_zdc; }

 protected:

  TowerInfo *_tmp_tower = nullptr;
  BbcPmtInfoV1 *_tmp_pmt = nullptr;
  CentralityInfo *_central = nullptr;
  BbcPmtInfoContainerV1 *_pmts_mbd = nullptr;
  TowerInfoContainer *_towers_zdc = nullptr;

  bool _zdc_check = false;
  unsigned int _key = std::numeric_limits<unsigned int>::max();

  int _op_mode = 0;
  int _use_ZDC = true;
  int _isMinBias = std::numeric_limits<int>::max();
  int _quality = std::numeric_limits<int>::max();
  int _side = std::numeric_limits<int>::max();
  int _channel = std::numeric_limits<int>::max();
  int _type = std::numeric_limits<int>::max();

  int _tubes_hit[2]{};
  int _tdc[2]{};

  float _mbd_charge_threshold = 0.4;

  float _zdc_energy_threshold = 1.;
  float _offset = std::numeric_limits<float>::signaling_NaN();
  float _energy = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum_n = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum_s = std::numeric_limits<float>::signaling_NaN();
  float _zdc_energy_sum = std::numeric_limits<float>::signaling_NaN();
  float _zdc_energy_sum_n = std::numeric_limits<float>::signaling_NaN();
  float _zdc_energy_sum_s = std::numeric_limits<float>::signaling_NaN();

  float _centrality_map[20]{};
  float _zdc_gain_factors[6]{};

  int m_mbd_side[128]{};
  int m_mbd_channel[128]{};
  float m_mbd_charge[128]{};
  float m_mbd_time_t[128]{};
  float m_mbd_time_q[128]{};
  float m_mbd_time[128]{};

  float m_zdc_energy_low[6]{};
  float m_zdc_energy_high[6]{};
  float m_zdc_sum_low[2]{};
  float m_zdc_sum_high[2]{};
};

#endif
