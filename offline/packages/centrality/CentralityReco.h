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

class BbcOutV1;


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
  BbcOutV1 *_bbc_out = nullptr;

  TowerInfoContainer *_towers_zdc = nullptr;
  CentralityInfo *_central = nullptr;

  bool _zdc_check = false;
  unsigned int _key = std::numeric_limits<unsigned int>::max();

  int _op_mode = 0;
  int _use_ZDC = true;
  int _isMinBias = std::numeric_limits<int>::max();
  int _quality = std::numeric_limits<int>::max();
  int _side = std::numeric_limits<int>::max();
  int _channel = std::numeric_limits<int>::max();
  int _type = std::numeric_limits<int>::max();

  float _mbd_charge_threshold = 0.4;

  short _tubes_hit_s;
  short _tubes_hit_n;

  float _zdc_energy_threshold = 1.;

  float _energy = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum_n = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum_s = std::numeric_limits<float>::signaling_NaN();
  float _zdc_energy_sum = std::numeric_limits<float>::signaling_NaN();
  float _zdc_energy_sum_n = std::numeric_limits<float>::signaling_NaN();
  float _zdc_energy_sum_s = std::numeric_limits<float>::signaling_NaN();

  float _centrality_map[20]{};
  float _zdc_gain_factors[6]{};

  float m_zdc_energy_low[6]{};
};

#endif
