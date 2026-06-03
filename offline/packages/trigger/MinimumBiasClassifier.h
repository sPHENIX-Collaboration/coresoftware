#ifndef TRIGGER_MINBIASCLASSIFIER_H
#define TRIGGER_MINBIASCLASSIFIER_H

#include "MinimumBiasInfo.h"

#include <fun4all/SubsysReco.h>

#include <array>
#include <limits>
#include <string>  // for allocator, string
#include <utility>
#include <vector>

// Forward declarations

class PHCompositeNode;
class Zdcinfo;
class MbdPmtContainer;
class MbdPmtHit;
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

  bool passesHitCut(MbdPmtHit *hit) const;

  int Download_centralityScale(const std::string &dbfile);
  int Download_centralityVertexScales(const std::string &dbfile);

  void setOverwriteScale(const std::string &url)
  {
    m_overwrite_url_scale = url;
  }
  void setOverwriteVtx(const std::string &url)
  {
    m_overwrite_url_vtx = url;
  }
  void setIsSim(const bool sim) { m_issim = sim; }

  void setSpecies(MinimumBiasInfo::SPECIES spec) { m_species = spec; };

  void abortEvents(const bool abort) { m_abortEvents = abort; };

  void set_minbiasNodeName(const std::string &name)
  {
    m_mb_info_nodename = name;
  }
  void set_mbdPmtNodeName(const std::string &name)
  {
    m_mbd_pmt_nodename = name;
  }
  void set_zdcInfoNodeName(const std::string &name)
  {
    m_zdc_info_nodename = name;
  }
  void set_globalvertexNodeName(const std::string &name)
  {
    m_global_vertex_nodename = name;
  }

 private:
  float getVertexScale();

  MinimumBiasInfo *m_mb_info{nullptr};
  MbdPmtContainer *m_mbd_container{nullptr};
  MbdPmtHit *m_mbd_pmt{nullptr};
  GlobalVertexMap *m_global_vertex_map{nullptr};
  Zdcinfo *m_zdcinfo{nullptr};

  bool m_abortEvents{false};
  bool m_issim{false};
  bool m_useZDC{true};
  bool m_box_cut{true};

  int m_hit_cut{2};

  double m_max_charge_cut{2100};
  double m_centrality_scale{std::numeric_limits<double>::quiet_NaN()};
  double m_vertex_scale{std::numeric_limits<double>::quiet_NaN()};

  float m_vertex{std::numeric_limits<float>::quiet_NaN()};

  static constexpr float m_z_vtx_cut{60.};
  static constexpr float m_mbd_north_cut{10.};
  static constexpr float m_mbd_south_cut{150};
  static constexpr float m_mbd_charge_cut{0.5};
  static constexpr float m_mbd_time_cut{25.};
  //  const int m_mbd_tube_cut{2};
  static constexpr float m_zdc_cut{60.};

  MinimumBiasInfo::SPECIES m_species{MinimumBiasInfo::SPECIES::AUAU};

  std::string m_mb_info_nodename{"MinimumBiasInfo"};
  std::string m_mbd_pmt_nodename{"MbdPmtContainer"};
  std::string m_zdc_info_nodename{"Zdcinfo"};
  std::string m_global_vertex_nodename{"GlobalVertexMap"};

  std::string m_overwrite_url_scale;
  std::string m_overwrite_url_vtx;


  std::array<float, 2> m_zdc_energy_sum{};
  std::array<float, 2> m_mbd_charge_sum{};
  std::array<int, 2> m_mbd_hit{};

  std::vector<std::pair<std::pair<float, float>, float>> m_vertex_scales{};
};

#endif
