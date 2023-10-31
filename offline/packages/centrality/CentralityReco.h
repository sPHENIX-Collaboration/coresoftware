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

class MbdOutV1;


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
  int FillCentralityInfo();
  int FillVars();
  int GetNodes(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);
  
  //! end of run method
  int End(PHCompositeNode *);

  void ResetVars();

 protected:

  MbdOutV1 *_mbd_out = nullptr;

  CentralityInfo *_central = nullptr;

  unsigned int _key = std::numeric_limits<unsigned int>::max();

  short _tubes_hit_s;
  short _tubes_hit_n;

  float _mbd_charge_sum = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum_n = std::numeric_limits<float>::signaling_NaN();
  float _mbd_charge_sum_s = std::numeric_limits<float>::signaling_NaN();

  float _centrality_map[20]{};
};

#endif
