#ifndef TRIGGER_MINBIASCLASSIFIER_H
#define TRIGGER_MINBIASCLASSIFIER_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <limits>
#include <string>  // for allocator, string
// Forward declarations

class MinimumBiasInfo;
class PHCompositeNode;
class MbdPmtContainer;
class MbdPmtHit;
class Zdcinfo;
class TowerInfo;
class GlobalVertexMap;
class GlobalVertex;

class MinimumBiasClassifier : public SubsysReco
{
 public:
  //! constructor
  explicit MinimumBiasClassifier(const std::string &name = "MinimumBiasClassifier");

  //! destructor

  ~MinimumBiasClassifier() override = default;

  int InitRun(PHCompositeNode *) override;
  void CreateNodes(PHCompositeNode *);
  int GetNodes(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *) override;
  //! end of run method

  int ResetEvent(PHCompositeNode *) override;

  int FillMinimumBiasInfo();

 private:
  const float _z_vtx_cut{60.};
  const float _mbd_north_cut{10.};
  const float _mbd_south_cut{150};
  const int _mbd_tube_cut{2};
  const float _zdc_cut{40.};
  const float charge_threshold = 0.4;
  const float time_threshold = 25.;

  MinimumBiasInfo *_mb_info{nullptr};
  Zdcinfo *_zdcinfo{nullptr};
  MbdPmtContainer *_mbd_pmt{nullptr};
  MbdPmtHit *_tmp_pmt{nullptr};
  GlobalVertexMap *_global_vertex_map{nullptr};

  std::array<float, 2> _zdc_energy_sum{};
  std::array<float, 2> _mbd_charge_sum{};
};

#endif
