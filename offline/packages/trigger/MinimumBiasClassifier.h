#ifndef TRIGGER_MINBIASCLASSIFIER_H
#define TRIGGER_MINBIASCLASSIFIER_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <limits>
#include <string>  // for allocator, string

// Forward declarations

class MinimumBiasInfo;
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
  ~MinimumBiasClassifier() override = default;

  //! full initialization
  int InitRun(PHCompositeNode *) override;
  void CreateNodes(PHCompositeNode *);
  int GetNodes(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *) override;
  int ResetEvent(PHCompositeNode *) override;

  int FillMinimumBiasInfo();

 private:

  const float _z_vtx_cut{60.};
  const float _mbd_north_cut{10.};
  const float _mbd_south_cut{150};
  const int _mbd_tube_cut{2};
  const float _zdc_cut{40.};

  MinimumBiasInfo *_mb_info{nullptr};
  MbdOut *_mbd_out{nullptr};
  GlobalVertexMap *_global_vertex_map{nullptr};
  TowerInfoContainer *_towers_zdc{nullptr};
  TowerInfo *_tmp_tower{nullptr};

  std::array<float, 2> _zdc_energy_sum{};
};

#endif
