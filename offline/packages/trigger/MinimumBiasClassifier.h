#ifndef MINBIASCLASSIFIER_H
#define MINBIASCLASSIFIER_H

#include "MinimumBiasInfov1.h"

#include <fun4all/SubsysReco.h>

#include <limits>

// Forward declarations

class MinimumBiasInfo;
class MinimumBiasInfov1;
class PHCompositeNode;
class MbdOut;
class TowerInfoContainer;
class TowerInfo;
class GlobalVertexMap;


class MinimumBiasClassifier : public SubsysReco
{
 public:
  //! constructor
  explicit MinimumBiasClassifier(const std::string &name = "MinimumBiasClassifier");

  //! destructor
  virtual ~MinimumBiasClassifier();

  //! full initialization
  int Init(PHCompositeNode *);
  int InitRun(PHCompositeNode *);
  void CreateNodes(PHCompositeNode *);
  int GetNodes(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);
  int FillMinimumBiasInfo();
  int FillVars();
  
  //! end of run method
  int End(PHCompositeNode *);

  virtual void Reset();


  
 private:
  const float _z_vtx_cut = 60.;
  const float _mbd_north_cut = 10.;
  const float _mbd_south_cut = 150;
  const int _mbd_tube_cut = 2;
  const float _zdc_cut = 40.;
  
  MinimumBiasInfov1 *_mb_info = nullptr;
  MbdOut *_mbd_out = nullptr;
  GlobalVertexMap *_global_vertex_map = nullptr;
  TowerInfoContainer *_towers_zdc = nullptr;
  TowerInfo *_tmp_tower = nullptr;
  float _energy = std::numeric_limits<float>::signaling_NaN();
  float _z_vertex = std::numeric_limits<float>::signaling_NaN();

  float _mbd_charge_sum[2] = {};

  float _zdc_energy_sum[2] = {};

  int _mbd_tubes_hit[2] = {};
};

#endif
