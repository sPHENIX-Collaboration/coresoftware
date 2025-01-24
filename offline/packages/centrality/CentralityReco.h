#ifndef CENTRALITY_CENTRALITYRECO_H
#define CENTRALITY_CENTRALITYRECO_H

#include <fun4all/SubsysReco.h>
#include <array>
#include <vector>
#include <limits>
#include <string>  // for string, allocator

// Forward declarations
class CentralityInfo;
class MinimumBiasInfo;
class PHCompositeNode;
class GlobalVertexMap;
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
    m_overwrite_divs = true;    
  }
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

  bool m_overwrite_divs{false};
  bool m_overwrite_scale{false};
  bool m_overwrite_vtx{false};
  std::string m_overwrite_url_divs{""};
  std::string m_overwrite_url_scale{""};
  std::string m_overwrite_url_vtx{""};

  const int NDIVS{100};
  const float mbd_charge_cut{0.5};
  const float mbd_time_cut{25};

  GlobalVertexMap *m_global_vertex_map{nullptr};
  MbdOut *m_mbd_out{nullptr};
  MbdPmtContainer *m_mbd_container{nullptr};
  MbdPmtHit *m_mbd_hit{nullptr};
  MinimumBiasInfo *m_mb_info{nullptr};
  CentralityInfo *m_central{nullptr};

  unsigned int m_key{std::numeric_limits<unsigned int>::max()};

  float m_mbd_total_charge{std::numeric_limits<float>::quiet_NaN()};

  double m_centrality_scale{std::numeric_limits<double>::quiet_NaN()};
  std::vector<std::pair<std::pair<float, float>, float>>  m_vertex_scales{};
  std::array<float, 100> m_centrality_map{};


};

#endif
