#ifndef CENTRALITY_CENTRALITYRECO_H
#define CENTRALITY_CENTRALITYRECO_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <limits>
#include <string>  // for string, allocator
#include <utility>
#include <vector>

// Forward declarations

class CentralityInfo;
class MinimumBiasInfo;
class PHCompositeNode;
class MbdOut;
class MbdPmtContainer;
class MbdPmtHit;

class CentralityReco : public SubsysReco
{
 public:
  //! constructor
  explicit CentralityReco(const std::string &name = "CentralityReco");

  //! destructor
  ~CentralityReco() override = default;

  //! full initialization
  int InitRun(PHCompositeNode *) override;
  void CreateNodes(PHCompositeNode *);
  int FillCentralityInfo();
  int FillVars();
  int GetNodes(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *) override;
  int ResetEvent(PHCompositeNode *) override;

  // Interface with CDB
  int Download_centralityDivisions(const std::string &dbfile);
  int Download_centralityScale(const std::string &dbfile);
  int Download_centralityVertexScales(const std::string &dbfile);

  void setOverwriteDivs(const std::string &url)
  {
    m_overwrite_url_divs = url;
  }
  void setOverwriteScale(const std::string &url)
  {
    m_overwrite_url_scale = url;
  }

  void setOverwriteVtx(const std::string &url)
  {
    m_overwrite_url_vtx = url;
  }

  void set_minbiasNodeName(const std::string &name)
  {
    m_mb_info_nodename = name;
  }
  void set_mbdOutNodeName(const std::string &name)
  {
    m_mbd_out_nodename = name;
  }
  void set_centralityNodeName(const std::string &name)
  {
    m_centrality_nodename = name;
  }
  void set_mbdPmtNodeName(const std::string &name)
  {
    m_mbd_pmt_nodename = name;
  }

 private:
  float getVertexScale();

  MbdOut *m_mbd_out{nullptr};
  MbdPmtContainer *m_mbd_container{nullptr};
  MbdPmtHit *m_mbd_hit{nullptr};
  MinimumBiasInfo *m_mb_info{nullptr};
  CentralityInfo *m_central{nullptr};

  static constexpr int NDIVS{100};

  static constexpr float mbd_charge_cut{0.5};
  static constexpr float mbd_time_cut{25};

  unsigned int m_key{std::numeric_limits<unsigned int>::max()};

  double m_centrality_scale{std::numeric_limits<double>::quiet_NaN()};

  float m_mbd_total_charge{0.};  // init to zero for use in first event


  std::string m_mb_info_nodename{"MinimumBiasInfo"};
  std::string m_mbd_out_nodename{"MbdOut"};
  std::string m_centrality_nodename{"CentralityInfo"};
  std::string m_mbd_pmt_nodename{"MbdPmtContainer"};

  std::string m_overwrite_url_divs;
  std::string m_overwrite_url_scale;
  std::string m_overwrite_url_vtx;

  std::vector<std::pair<std::pair<float, float>, float>> m_vertex_scales{};
  std::array<float, 100> m_centrality_map{};
};

#endif
