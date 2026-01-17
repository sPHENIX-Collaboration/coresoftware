#include "StateClusterResidualsQA.h"

#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/ActsGeometry.h>

#include <trackbase_historic/SvtxAlignmentState.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <Rtypes.h>
#include <TH1F.h>

#include <format>

namespace
{
  template <typename T>
  inline T square (T const& t) { return t * t; }

  template <class T>
  class range_adaptor
  {
   public:
    explicit range_adaptor(
        T const& begin,
        T const& end)
      : m_begin(begin)
      , m_end(end)
    {
    }
    T const& begin() { return m_begin; }
    T const& end() { return m_end; }

   private:
    T m_begin;
    T m_end;
  };
}  // namespace

StateClusterResidualsQA::StateClusterResidualsQA(const std::string& name)
  : SubsysReco(name)
{
}

int StateClusterResidualsQA::InitRun(
    PHCompositeNode* top_node)
{
  createHistos();

  // F4A will not actually ABORTRUN unless that return code is issued here
  auto* track_map = findNode::getClass<SvtxTrackMap>(top_node, m_track_map_node_name);
  if (!track_map)
  {
    std::cout
        << PHWHERE << "\n"
        << "\tCould not get track map:\n"
        << "\t\"" << m_track_map_node_name << "\"\n"
        << "\tAborting\n"
        << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto* cluster_map = findNode::getClass<TrkrClusterContainer>(top_node, m_clusterContainerName);
  if (!cluster_map)
  {
    std::cout
        << PHWHERE << "\n"
        << "\tCould not get cluster map:\n"
        << "\t\"" << m_clusterContainerName << "\"\n"
        << "\tAborting\n"
        << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  auto *geometry = findNode::getClass<ActsGeometry>(top_node, "ActsGeometry");
  if (!geometry)
  {
    std::cout
        << PHWHERE << "\n"
        << "\tCould not get ActsGeometry:\n"
        << "\t\"" << "ActsGeometry" << "\"\n"
        << "\tAborting\n"
        << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto* hm = QAHistManagerDef::getHistoManager();
  if (!hm)
  {
    std::cout
        << PHWHERE << "\n"
        << "\tCould not get QAHistManager\n"
        << "\tAborting\n"
        << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (const auto& cfg : m_pending)
  {
    m_histograms_x.push_back(dynamic_cast<TH1 *>(hm->getHisto(std::string(cfg.name + "_x"))));
    m_histograms_y.push_back(dynamic_cast<TH1 *>(hm->getHisto(std::string(cfg.name + "_y"))));
    m_histograms_z.push_back(dynamic_cast<TH1 *>(hm->getHisto(std::string(cfg.name + "_z"))));
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int StateClusterResidualsQA::process_event(PHCompositeNode* top_node)
{
  auto* track_map = findNode::getClass<SvtxTrackMap>(top_node, m_track_map_node_name);
  auto *cluster_map = findNode::getClass<TrkrClusterContainer>(top_node, m_clusterContainerName);
  auto *geometry = findNode::getClass<ActsGeometry>(top_node, "ActsGeometry");

  for (auto const& [idkey, track] : *track_map)
  {
    if (!track)
    {
      continue;
    }

    // count states
    std::map<TrkrDefs::TrkrId, int> counters = {
        {TrkrDefs::mvtxId, 0},
        {TrkrDefs::inttId, 0},
        {TrkrDefs::tpcId, 0},
        {TrkrDefs::micromegasId, 0},
    };

    for (auto const& [path_length, state] : range_adaptor(track->begin_states(), track->end_states()))
    {
      // There is an additional state representing the vertex at the beginning of the map,
      // but getTrkrId will return 0 for its corresponding cluster
      // Identify it as having path_length identically equal to 0
      if (path_length == 0) { continue; }

      auto trkr_id = static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(state->get_cluskey()));
      auto itr = counters.find(trkr_id);
      if (itr == counters.end()) { continue; }
      ++itr->second;
    }

    float track_eta = track->get_eta();
    float track_phi = track->get_phi();
    float track_pt = track->get_pt();
    int h = 0;
    for (const auto& cfg : m_pending)
    {
      if (cfg.charge != 0)
      {
        if ((cfg.charge < 0) && track->get_positive_charge())
        {
          continue;
        }
        else if ((cfg.charge > 0) && !(track->get_positive_charge()))
        {
          continue;
        }
      }
      if (cfg.min_mvtx_clusters <= counters[TrkrDefs::mvtxId] && cfg.max_mvtx_clusters >= counters[TrkrDefs::mvtxId]
          && cfg.min_intt_clusters <= counters[TrkrDefs::inttId] && cfg.max_intt_clusters >= counters[TrkrDefs::inttId]
          && cfg.min_tpc_clusters <= counters[TrkrDefs::tpcId] && cfg.max_tpc_clusters >= counters[TrkrDefs::tpcId]
          && cfg.phi_min <= track_phi && cfg.phi_max >= track_phi
          && cfg.eta_min <= track_eta && cfg.eta_max >= track_eta
          && cfg.pt_min <= track_pt && cfg.pt_max >= track_pt)
      {
        for (auto const& [path_length, state] : range_adaptor(track->begin_states(), track->end_states()))
        {
          if (path_length == 0) { continue; }

          auto *cluster = cluster_map->findCluster(state->get_cluskey());
          float state_x = state->get_x();
          float state_y = state->get_y();
          float state_z = state->get_z();
          Acts::Vector3 glob = geometry->getGlobalPosition(state->get_cluskey(), cluster);
          float cluster_x = glob.x();
          float cluster_y = glob.y();
          float cluster_z = glob.z();
          if (cluster)
          {
            m_histograms_x[h]->Fill(state_x - cluster_x);
            m_histograms_y[h]->Fill(state_y - cluster_y);
            m_histograms_z[h]->Fill(state_z - cluster_z);
          }
        }
      } 
      ++h;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void StateClusterResidualsQA::createHistos()
{
  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for (const auto& cfg : m_pending)
  {
    TH1F* h_new_x = new TH1F(
        (cfg.name + "_x").c_str(),
        ";State-Cluster X Residual [cm];Entries",
        m_nBins, m_xrange.first, m_xrange.second);
    h_new_x->SetMarkerColor(kBlue);
    h_new_x->SetLineColor(kBlue);
    hm->registerHisto(h_new_x);
    TH1F* h_new_y = new TH1F(
        (cfg.name + "_y").c_str(),
        ";State-Cluster Y Residual [cm];Entries",
        m_nBins, m_yrange.first, m_yrange.second);
    h_new_y->SetMarkerColor(kBlue);
    h_new_y->SetLineColor(kBlue);
    hm->registerHisto(h_new_y);
    TH1F* h_new_z = new TH1F(
        (cfg.name + "_z").c_str(),
        ";State-Cluster Z Residual [cm];Entries",
        m_nBins, m_zrange.first, m_zrange.second);
    h_new_z->SetMarkerColor(kBlue);
    h_new_z->SetLineColor(kBlue);
    hm->registerHisto(h_new_z);
  }
}

int StateClusterResidualsQA::EndRun(const int /*unused*/)
{
  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  return Fun4AllReturnCodes::EVENT_OK;
}
