
#include "TrackQA.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <qautils/QAHistManagerDef.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <TH2.h>
#include <trackbase/ActsGeometry.h>
//____________________________________________________________________________..
TrackQA::TrackQA(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TrackQA::InitRun(PHCompositeNode* /*unused*/)
{
  createHistos();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackQA::process_event(PHCompositeNode* topNode)
{
  auto clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  auto geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  auto trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!trackmap or !clustermap or !geometry)
  {
    std::cout << PHWHERE << "Missing node(s), can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  auto h_ntrack = dynamic_cast<TH2*>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks").c_str()));
  auto h_nmaps = dynamic_cast<TH1*>(hm->getHisto(std::string(getHistoPrefix() + "nmaps").c_str()));
  auto h_nintt = dynamic_cast<TH1*>(hm->getHisto(std::string(getHistoPrefix() + "nintt").c_str()));
  auto h_ntpc = dynamic_cast<TH1*>(hm->getHisto(std::string(getHistoPrefix() + "ntpc").c_str()));
  auto h_nmms = dynamic_cast<TH1*>(hm->getHisto(std::string(getHistoPrefix() + "ntpot").c_str()));
  auto h_lxresid = dynamic_cast<TH2*>(hm->getHisto(std::string(getHistoPrefix() + "lxlineresiduals").c_str()));
  auto h_lzresid = dynamic_cast<TH2*>(hm->getHisto(std::string(getHistoPrefix() + "lzlineresiduals").c_str()));
  auto h_gzresid = dynamic_cast<TH2*>(hm->getHisto(std::string(getHistoPrefix() + "gzlineresiduals").c_str()));
  auto h_gyresid = dynamic_cast<TH2*>(hm->getHisto(std::string(getHistoPrefix() + "gylineresiduals").c_str()));
  auto h_gxresid = dynamic_cast<TH2*>(hm->getHisto(std::string(getHistoPrefix() + "gxlineresiduals").c_str()));

  m_tracks += trackmap->size();
  for (const auto& [key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }

    auto ckeys = get_cluster_keys(track);
    std::vector<Acts::Vector3> cluspos;
    TrackFitUtils::getTrackletClusters(geometry, clustermap, cluspos, ckeys);

    auto lineFitParams = lineFitClusters(cluspos);

    /// Only look at tracks that go through 20 cm of origin
    if (std::fabs(std::get<1>(lineFitParams)) > 20)
    {
      continue;
    }
    float px = track->get_px();
    float py = track->get_py();
    float pz = track->get_pz();
    float pt = std::sqrt(QAG4Util::square(px) + QAG4Util::square(py));
    float eta = std::atanh(pz / std::sqrt(QAG4Util::square(pt) + QAG4Util::square(pz)));
    float phi = std::atan2(py, px);
    h_ntrack->Fill(phi, eta);

    int nmaps = 0;
    int nintt = 0;
    int ntpc = 0;
    int nmms = 0;

    for (auto& ckey : ckeys)
    {
      switch (TrkrDefs::getTrkrId(ckey))
      {
      case TrkrDefs::mvtxId:
        nmaps++;
        break;
      case TrkrDefs::inttId:
        nintt++;
        break;
      case TrkrDefs::tpcId:
        ntpc++;
        break;
      case TrkrDefs::micromegasId:
        nmms++;
        break;
      }
    }
    if (nmaps == 0)
    {
      continue;
    }
    int i = 0;
    for (auto& ckey : ckeys)
    {
      auto& glob = cluspos[i];
      i++;
      auto cluster = clustermap->findCluster(ckey);

      auto intersection = TrackFitUtils::surface_3Dline_intersection(ckey, cluster, geometry,
                                                                     std::get<0>(lineFitParams), std::get<1>(lineFitParams), std::get<2>(lineFitParams), std::get<3>(lineFitParams));

      auto surf = geometry->maps().getSurface(ckey, cluster);

      Acts::Vector3 surfnorm = surf->normal(geometry->geometry().getGeoContext());
      float statelx = std::numeric_limits<float>::quiet_NaN();
      float statelz = std::numeric_limits<float>::quiet_NaN();
      if (!std::isnan(intersection.x()))
      {
        auto locstateres = surf->globalToLocal(geometry->geometry().getGeoContext(),
                                               intersection * Acts::UnitConstants::cm,
                                               surfnorm);
        if (locstateres.ok())
        {
          Acts::Vector2 loc = locstateres.value() / Acts::UnitConstants::cm;
          statelx = loc(0);
          statelz = loc(1);
        }
        else
        {
          Acts::Vector3 loc = surf->transform(geometry->geometry().getGeoContext()).inverse() * (intersection * Acts::UnitConstants::cm);
          loc /= Acts::UnitConstants::cm;
          statelx = loc(0);
          statelz = loc(1);
        }
      }
      unsigned int layer = TrkrDefs::getLayer(ckey);
      auto localclus = geometry->getLocalCoords(ckey, cluster);

      h_lxresid->Fill(layer, localclus.x() - statelx);
      h_lzresid->Fill(layer, localclus.y() - statelz);
      h_gxresid->Fill(layer, glob.x() - intersection.x());
      h_gyresid->Fill(layer, glob.y() - intersection.y());
      h_gzresid->Fill(layer, glob.z() - intersection.z());
    }
    h_nmaps->Fill(nmaps);
    h_nintt->Fill(nintt);
    h_ntpc->Fill(ntpc);
    h_nmms->Fill(nmms);
  }

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}
std::vector<TrkrDefs::cluskey> TrackQA::get_cluster_keys(SvtxTrack* track)
{
  std::vector<TrkrDefs::cluskey> out;
  for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
  {
    if (seed)
    {
      std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
    }
  }
  return out;
}

std::tuple<float, float, float, float>
TrackQA::lineFitClusters(std::vector<Acts::Vector3>& positions) const
{
  TrackFitUtils::position_vector_t xypoints, rzpoints;
  for (auto& pos : positions)
  {
    float clusr = std::sqrt(QAG4Util::square(pos.x()) + QAG4Util::square(pos.y()));
    if (pos.y() < 0)
    {
      clusr *= -1;
    }
    // exclude silicon and tpot clusters for now
    if (std::fabs(clusr) > 80 || std::fabs(clusr) < 30)
    {
      continue;
    }
    rzpoints.push_back(std::make_pair(pos.z(), clusr));
    xypoints.push_back(std::make_pair(pos.x(), pos.y()));
  }

  auto xyparams = TrackFitUtils::line_fit(xypoints);
  auto rzparams = TrackFitUtils::line_fit(rzpoints);

  return std::make_tuple(std::get<0>(xyparams),
                         std::get<1>(xyparams),
                         std::get<0>(rzparams),
                         std::get<1>(rzparams));
}
//____________________________________________________________________________..
int TrackQA::EndRun(const int runnumber)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH2* h_tracksperevent = dynamic_cast<TH2*>(hm->getHisto(std::string(getHistoPrefix() + "ntracksperrun").c_str()));
  h_tracksperevent->Fill(runnumber, m_tracks / m_event);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackQA::End(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string TrackQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void TrackQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpc").c_str(),
                      "TPC clusters per track", 150, 0, 150);
    h->GetXaxis()->SetTitle("nTPC");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nmaps").c_str(),
                      "MVTX clusters per track", 20, 0, 20);
    h->GetXaxis()->SetTitle("nMVTX");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nintt").c_str(),
                      "INTT clusters per track", 20, 0, 20);
    h->GetXaxis()->SetTitle("nINTT");
    hm->registerHisto(h);
  }
  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpot").c_str(),
                      "TPOT clusters per track", 6, 0, 6);
    h->GetXaxis()->SetTitle("nTPOT");
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "lxlineresiduals").c_str(),
                      "Local x residuals", 57, 0, 57, 1000, -10, 10);
    h->GetXaxis()->SetTitle("Layer");
    h->GetYaxis()->SetTitle("l_{x}^{clus}-l_{x}^{state} [cm]");
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "lzlineresiduals").c_str(),
                      "Local z residuals", 57, 0, 57, 1000, -10, 10);
    h->GetXaxis()->SetTitle("Layer");
    h->GetYaxis()->SetTitle("l_{z}^{clus}-l_{z}^{state} [cm]");
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "gzlineresiduals").c_str(),
                      "Global z residuals", 57, 0, 57, 1000, -10, 10);
    h->GetXaxis()->SetTitle("Layer");
    h->GetYaxis()->SetTitle("g_{z}^{clus}-g_{z}^{state} [cm]");
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "gylineresiduals").c_str(),
                      "Global y residuals", 57, 0, 57, 1000, -10, 10);
    h->GetXaxis()->SetTitle("Layer");
    h->GetYaxis()->SetTitle("g_{y}^{clus}-g_{y}^{state} [cm]");
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "gxlineresiduals").c_str(),
                      "Global x residuals", 57, 0, 57, 1000, -10, 10);
    h->GetXaxis()->SetTitle("Layer");
    h->GetYaxis()->SetTitle("g_{x}^{clus}-g_{x}^{state} [cm]");
    hm->registerHisto(h);
  }
  {
    auto h = new TH2F(std::string(getHistoPrefix() + "nrecotracks").c_str(),
                      "Num reconstructed tracks", 300, -3.14159, 3.1459, 100, -1.1, 1.1);
    h->GetXaxis()->SetTitle("#phi [rad]");
    h->GetYaxis()->SetTitle("#eta");
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "ntracksperrun").c_str(),
                      "Num reconstructed tracks per run", m_runbins, m_beginRun, m_endRun, 100, 0, 1);
    h->GetYaxis()->SetTitle("N_{tracks}/event");
    h->GetXaxis()->SetTitle("Run number");
    hm->registerHisto(h);
  }
}