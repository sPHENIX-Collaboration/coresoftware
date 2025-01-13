#ifndef TRIGGER_MINBIASCLASSIFIER_H
#define TRIGGER_MINBIASCLASSIFIER_H

#include <fun4all/SubsysReco.h>
#include <vector>
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
class MbdPmtContainer;
class MbdPmtHit;
class Zdcinfo;
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

  bool passesHitCut(MbdPmtHit *hit);

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

 private:
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
  std::vector<std::pair<std::pair<float, float>, float>>  m_vertex_scales{};

};

#endif
