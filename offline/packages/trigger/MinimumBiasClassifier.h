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

  //! bunch class of a t0 relative to the run's main-bunch offset (MinimumBiasInfo::BUNCH)
  int bunchClass(float t0_corrected) const;

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

  void setRejectPileup(bool v) { m_reject_pileup = v; };

  //! Au+Au data: bunch classes (main, -5 ns, +5 ns) are used when the CentralityScale payload carries a t0
  //! offset; false ignores it and keeps the single-class selection
  void setUseBunchClasses(bool v) { m_use_bunch_classes = v; }
  bool bunchClassesActive() const { return m_bunch_classes; }

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
  float getVertexScale(int bunch) const;

  MinimumBiasInfo *m_mb_info{nullptr};
  MbdPmtContainer *m_mbd_container{nullptr};
  MbdPmtHit *m_mbd_pmt{nullptr};
  GlobalVertexMap *m_global_vertex_map{nullptr};
  Zdcinfo *m_zdcinfo{nullptr};

  bool m_abortEvents{false};
  bool m_issim{false};
  bool m_useZDC{true};
  bool m_box_cut{true};
  bool m_reject_pileup{true};
  bool m_use_bunch_classes{true};
  bool m_bunch_classes{false};  // this run: bunch-class payload found and used

  float m_pileup_charge_cut{200.};

  int m_hit_cut{2};

  double m_max_charge_cut{2100};
  // per bunch class (MAIN, M5, P5); single-class payloads fill only MAIN
  std::array<double, 3> m_centrality_scale{};
  double m_vertex_scale{std::numeric_limits<double>::quiet_NaN()};

  float m_vertex{std::numeric_limits<float>::quiet_NaN()};

  float m_z_vtx_cut{60.};
  float m_mbd_north_cut{10.};
  float m_mbd_south_cut{150.};
  float m_mbd_charge_cut{0.5};
  float m_mbd_time_cut{25.};
  //  const int m_mbd_tube_cut{2};
  float m_zdc_cut{60.};

  // bunch-class selection (Au+Au data with a bunch-class payload), matching the calibration
  float m_bunch_zdc_cut{40.};
  float m_bunch_charge_cut{0.4};
  double m_t0_offset{0.};
  double m_spacing_m{std::numeric_limits<double>::quiet_NaN()};
  double m_spacing_p{std::numeric_limits<double>::quiet_NaN()};
  // per-tube time cut, |t_hit - t0 offset| < this (the bunch windows themselves are applied to the event t0)
  double m_hit_time_window{25.};

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

  std::array<std::vector<std::pair<std::pair<float, float>, float>>, 3> m_vertex_scales{};
};

#endif
