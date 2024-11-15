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
class MbdPmtContainer;
class MbdPmtHit;
class Zdcinfo;
class GlobalVertexMap;

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

  bool passesHitCut(MbdPmtHit *hit);

 private:
  const float m_z_vtx_cut{60.};
  const float m_mbd_north_cut{10.};
  const float m_mbd_south_cut{150};
  const float m_mbd_charge_cut{0.4};
  const float m_mbd_time_cut{20.};
//  const int m_mbd_tube_cut{2};
  const float m_zdc_cut{40.};

  MinimumBiasInfo *m_mb_info{nullptr};
  MbdOut *m_mbd_out{nullptr};
  MbdPmtContainer *m_mbd_container{nullptr};
  MbdPmtHit *m_mbd_pmt{nullptr};
  GlobalVertexMap *m_global_vertex_map{nullptr};
  Zdcinfo *m_zdcinfo{nullptr};

  std::array<float, 2> m_zdc_energy_sum{};
  std::array<float, 2> m_mbd_charge_sum{};
  std::array<int, 2> m_mbd_hit{};
};

#endif
