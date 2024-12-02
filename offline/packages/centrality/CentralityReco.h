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

 private:


  float getVertexScale();

  std::string m_dbfilename;

  const int NDIVS{100};
  const float mbd_charge_cut{0.5};
  const float mbd_time_cut{20};

  GlobalVertexMap *m_global_vertex_map{nullptr};
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
