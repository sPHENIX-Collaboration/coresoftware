#include "TrackFittingQA.h"

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
#include <TH2F.h>
#include <TProfile.h>

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

TrackFittingQA::TrackFittingQA(const std::string& name)
  : SubsysReco(name)
{
  // ...
}

int TrackFittingQA::Init(PHCompositeNode* /*unused*/)
{
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

  // Create histograms

  for (int charge = 0; charge < 2; ++charge)
  {
    std::string charge_str = charge ? "positively" : "negatively";
    EColor color = charge ? kRed : kBlue;

    delete m_quality_hist[charge];
    m_quality_hist[charge] = new TH1F(
        std::format("h_{}_quality_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("quality distribution ({} charged tracks);quality;Counts", charge_str).c_str(),
        50, m_quality_hist_xrange.first, m_quality_hist_xrange.second);
    hm->registerHisto(m_quality_hist[charge]);
    m_quality_hist[charge]->SetMarkerColor(color);
    m_quality_hist[charge]->SetLineColor(color);
    // ...

    delete m_p_hist[charge];
    m_p_hist[charge] = new TH1F(
        std::format("h_{}_p_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("p distribution ({} charged tracks);p (GeV);Counts", charge_str).c_str(),
        50, m_p_hist_xrange.first, m_p_hist_xrange.second);
    hm->registerHisto(m_p_hist[charge]);
    m_p_hist[charge]->SetMarkerColor(color);
    m_p_hist[charge]->SetLineColor(color);
    // ...

    delete m_pt_hist[charge];
    m_pt_hist[charge] = new TH1F(
        std::format("h_{}_pt_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("pt distribution ({} charged tracks);pt (GeV);Counts", charge_str).c_str(),
        50, m_pt_hist_xrange.first, m_pt_hist_xrange.second);
    hm->registerHisto(m_pt_hist[charge]);
    m_pt_hist[charge]->SetMarkerColor(color);
    m_pt_hist[charge]->SetLineColor(color);
    // ...

    delete m_pt_err_hist[charge];
    m_pt_err_hist[charge] = new TH2F(
        std::format("h_{}_pt_err_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("relative pt err vs pt distribution ({} charged tracks);pt (GeV);pt err / pt (unitless)", charge_str).c_str(),
        50, m_pt_err_hist_xrange.first, m_pt_err_hist_xrange.second,
        50, m_pt_err_hist_yrange.first, m_pt_err_hist_yrange.second);
    hm->registerHisto(m_pt_err_hist[charge]);
    m_pt_err_hist[charge]->SetMarkerColor(color);
    m_pt_err_hist[charge]->SetLineColor(color);
    // ...

    delete m_pt_err_profile[charge];
    m_pt_err_profile[charge] = new TProfile(
        std::format("h_{}_pt_err_{}_charged_tracks_profile", Name(), charge_str).c_str(),
        std::format("relative pt err vs pt distribution ({} charged tracks);pt (GeV);pt err / pt (unitless)", charge_str).c_str(),
        50, m_pt_err_hist_xrange.first, m_pt_err_hist_xrange.second,
        m_pt_err_hist_yrange.first, m_pt_err_hist_yrange.second);
    hm->registerHisto(m_pt_err_profile[charge]);
    m_pt_err_profile[charge]->SetMarkerColor(color);
    m_pt_err_profile[charge]->SetLineColor(color);
    // ...

    delete m_eta_hist[charge];
    m_eta_hist[charge] = new TH1F(
        std::format("h_{}_eta_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("eta distribution ({} charged tracks);eta;Counts", charge_str).c_str(),
        50, -1.2, 1.2);
    hm->registerHisto(m_eta_hist[charge]);
    m_eta_hist[charge]->SetMarkerColor(color);
    m_eta_hist[charge]->SetLineColor(color);
    // ...

    delete m_phi_eta_hist[charge];
    m_phi_eta_hist[charge] = new TH2F(
        std::format("h_{}_phi_eta_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("phi-eta distribution ({} charged tracks);phi (radians);eta", charge_str).c_str(),
        50, -3.2, 3.2,
        50, -1.2, 1.2);
    hm->registerHisto(m_phi_eta_hist[charge]);
    // ...

    delete m_mvtx_states_hist[charge];
    m_mvtx_states_hist[charge] = new TH1F(
        std::format("h_{}_mvtx_states_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("mvtx states distribution ({} charged tracks);number of mvtx states;tracks with that number of mvtx states", charge_str).c_str(),
        6, -0.5, 5.5);
    hm->registerHisto(m_mvtx_states_hist[charge]);
    m_mvtx_states_hist[charge]->SetMarkerColor(color);
    m_mvtx_states_hist[charge]->SetLineColor(color);
    // ...

    delete m_intt_states_hist[charge];
    m_intt_states_hist[charge] = new TH1F(
        std::format("h_{}_intt_states_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("intt states distribution ({} charged tracks);number of intt states;tracks with that number of intt states", charge_str).c_str(),
        5, -0.5, 4.5);
    hm->registerHisto(m_intt_states_hist[charge]);
    m_intt_states_hist[charge]->SetMarkerColor(color);
    m_intt_states_hist[charge]->SetLineColor(color);
    // ...

    delete m_tpc_states_hist[charge];
    m_tpc_states_hist[charge] = new TH1F(
        std::format("h_{}_tpc_states_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("tpc states distribution ({} charged tracks);number of tpc states;tracks with that number of tpc states", charge_str).c_str(),
        61, -0.5, 60.5);
    hm->registerHisto(m_tpc_states_hist[charge]);
    m_tpc_states_hist[charge]->SetMarkerColor(color);
    m_tpc_states_hist[charge]->SetLineColor(color);
    // ...

    delete m_tpot_states_hist[charge];
    m_tpot_states_hist[charge] = new TH1F(
        std::format("h_{}_tpot_states_{}_charged_tracks", Name(), charge_str).c_str(),
        std::format("tpot states distribution ({} charged tracks);number of tpot states;tracks with that number of tpot states", charge_str).c_str(),
        4, -0.5, 3.5);
    hm->registerHisto(m_tpot_states_hist[charge]);
    m_tpot_states_hist[charge]->SetMarkerColor(color);
    m_tpot_states_hist[charge]->SetLineColor(color);
    // ...
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackFittingQA::InitRun(
    PHCompositeNode* top_node)
{
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

  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackFittingQA::process_event(
    PHCompositeNode* top_node)
{
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

  for (auto const& [idkey, track] : *track_map)
  {
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

    // Cuts
    if ( track->get_quality() < m_min_quality ) { continue; }
    if ( m_max_quality < track->get_quality() ) { continue; }
    if ( track->get_p() < m_min_p ) { continue; }
    if ( track->get_pt() < m_min_pt ) { continue; }
    if ( m_max_abs_eta < std::abs(track->get_eta()) ) { continue; }
    if ( counters[TrkrDefs::inttId] < m_min_intt_states ) { continue; }
    if ( counters[TrkrDefs::mvtxId] < m_min_mvtx_states ) { continue; }
    if ( counters[TrkrDefs::tpcId] < m_min_tpc_states ) { continue; }
    if ( counters[TrkrDefs::micromegasId] < m_min_tpot_states ) { continue; }
    if ( track->get_crossing() < m_min_crossing ) { continue; }
    if ( m_max_crossing < track->get_crossing() ) { continue; }

    // Fill histograms if passing
    int charge = (0 < track->get_charge()) ? 1 : 0;

    // Immitated from TrackResiduals.cc:1896-1899 (as of commit b75c9ef557d096f7d6ec9bbf2082d4291eddc4a1)
    float px = track->get_px();
    float py = track->get_py();
    float CVxx = track->get_error(3, 3);
    float CVxy = track->get_error(3, 4);
    float CVyy = track->get_error(4, 4);
    float pt_err = std::sqrt((CVxx * square(px) + 2 * CVxy * px * py + CVyy * square(py)) / (square(px) + square(py)));

    m_quality_hist[charge]->Fill(track->get_quality());
    m_p_hist[charge]->Fill(track->get_p());
    m_pt_hist[charge]->Fill(track->get_pt());
    m_pt_err_hist[charge]->Fill(track->get_pt(), pt_err);
    m_pt_err_profile[charge]->Fill(track->get_pt(), pt_err);
    m_eta_hist[charge]->Fill(track->get_eta());
    m_phi_eta_hist[charge]->Fill(track->get_phi(), track->get_eta());

    m_mvtx_states_hist[charge]->Fill(counters[TrkrDefs::mvtxId]);
    m_intt_states_hist[charge]->Fill(counters[TrkrDefs::inttId]);
    m_tpc_states_hist[charge]->Fill(counters[TrkrDefs::tpcId]);
    m_tpot_states_hist[charge]->Fill(counters[TrkrDefs::micromegasId]);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackFittingQA::ResetEvent(
    PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackFittingQA::EndRun(
    int const /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackFittingQA::End(
    PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrackFittingQA::Reset(
    PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrackFittingQA::Print(
    std::string const& /*unused*/) const
{
  // ...
}
