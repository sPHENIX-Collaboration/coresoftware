#ifndef TRIGGER_MINBIASCLASSIFIER_H
#define TRIGGER_MINBIASCLASSIFIER_H

#include <phparameter/PHParameters.h>
#include <fun4all/SubsysReco.h>
#include <array>
#include <limits>
#include <string>  // for allocator, string
#include <utility>
#include <vector>

#include "MinimumBiasInfo.h"
// Forward declarations

class MinimumBiasInfo;
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
  static void CreateNodes(PHCompositeNode *);
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
    m_overwrite_scale = true;
  }
  void setOverwriteVtx(const std::string &url)
  {
    m_overwrite_url_vtx = url;
    m_overwrite_vtx = true;
  }
  void setIsSim(const bool sim) { m_issim = sim; }
  void set_mbd_total_charge_cut(const double max_charge_cut) { m_max_charge_cut = max_charge_cut; }

  void setSpecies(MinimumBiasInfo::SPECIES spec) { m_species = spec; };
  
 private:
  bool m_issim{false};
  bool m_useZDC{true};
  bool m_box_cut{true};
  int m_hit_cut{2};
  double m_max_charge_cut{2100};
  
  MinimumBiasInfo::SPECIES m_species{MinimumBiasInfo::SPECIES::AUAU};

  float getVertexScale();
  std::string m_dbfilename;

  bool m_overwrite_scale{false};
  bool m_overwrite_vtx{false};
  std::string m_overwrite_url_scale{""};
  std::string m_overwrite_url_vtx{""};

  const float m_z_vtx_cut{60.};
  const float m_mbd_north_cut{10.};
  const float m_mbd_south_cut{150};
  const float m_mbd_charge_cut{0.5};
  const float m_mbd_time_cut{25.};
  //  const int m_mbd_tube_cut{2};
  const float m_zdc_cut{60.};

  MinimumBiasInfo *m_mb_info{nullptr};
  MbdPmtContainer *m_mbd_container{nullptr};
  MbdPmtHit *m_mbd_pmt{nullptr};
  GlobalVertexMap *m_global_vertex_map{nullptr};
  Zdcinfo *m_zdcinfo{nullptr};

  std::array<float, 2> m_zdc_energy_sum{};
  std::array<float, 2> m_mbd_charge_sum{};
  std::array<int, 2> m_mbd_hit{};

  double m_centrality_scale{std::numeric_limits<double>::quiet_NaN()};
  double m_vertex_scale{std::numeric_limits<double>::quiet_NaN()};
  float m_vertex{std::numeric_limits<float>::quiet_NaN()};
  std::vector<std::pair<std::pair<float, float>, float>> m_vertex_scales{};

  PHParameters m_MinBiasParams;
  PHCompositeNode* m_parNode{nullptr};
};

#endif
