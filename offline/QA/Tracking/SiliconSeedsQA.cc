
#include "SiliconSeedsQA.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <qautils/QAHistManagerDef.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrackFitUtils.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include <trackbase/ActsGeometry.h>
//____________________________________________________________________________..
SiliconSeedsQA::SiliconSeedsQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int SiliconSeedsQA::InitRun(PHCompositeNode * /*unused*/)
{
  createHistos();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SiliconSeedsQA::process_event(PHCompositeNode *topNode)
{
  auto clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  auto geometry = findNode::getClass<ActsGeometry>(topNode, m_actsgeometryName);
  auto trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);

  if (!trackmap or !clustermap or !geometry or !vertexmap)
  {
    std::cout << PHWHERE << "Missing node(s), can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  h_ntrack1d->Fill(trackmap->size());

  std::pair<int, int> ntrack_isfromvtx;    // first: number of tracks associated to a vertex, second: number of tracks not associated to a vertex
  std::pair<int, int> ntrack_isposcharge;  // first: number of tracks with negative charge, second: number of tracks with positive charge

  for (const auto &[key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }

    auto ckeys = get_cluster_keys(track);
    float eta = track->get_eta();
    float phi = track->get_phi();
    float pt = track->get_pt();
    int charge = track->get_charge();
    float chi2ndf = track->get_quality();

    int trkcrossing = track->get_crossing();
    // std::cout << "track crossing: " << trkcrossing << std::endl;

    int nmaps = 0;
    int nintt = 0;
//    int ntpc = 0;
//    int nmms = 0;

    for (auto &ckey : ckeys)
    {
      switch (TrkrDefs::getTrkrId(ckey))
      {
      case TrkrDefs::mvtxId:
        nmaps++;
        break;
      case TrkrDefs::inttId:
        nintt++;
        break;
      // case TrkrDefs::tpcId:
      //   ntpc++;
      //   break;
      // case TrkrDefs::micromegasId:
      //   nmms++;
      //   break;
      default:
	break;
      }
    }

    Acts::Vector3 zero = Acts::Vector3::Zero();
    auto dcapair_origin = TrackAnalysisUtils::get_dca(track, zero);

    auto trackvtx = vertexmap->get(track->get_vertex_id());
    if (!trackvtx)
    {
      ntrack_isfromvtx.first++;
      continue;
    }
    ntrack_isfromvtx.second++;

    if (charge > 0)
    {
      ntrack_isposcharge.second++;
    }
    else
    {
      ntrack_isposcharge.first++;
    }

    Acts::Vector3 track_vtx(trackvtx->get_x(), trackvtx->get_y(), trackvtx->get_z());
    auto dcapair_vtx = TrackAnalysisUtils::get_dca(track, track_vtx);

    h_ntrack->Fill(eta, phi);
    h_nmaps->Fill(nmaps);
    h_nintt->Fill(nintt);
    h_nmaps_nintt->Fill(nmaps, nintt);
    h_avgnclus_eta_phi->Fill(eta, phi, nmaps + nintt);
    h_trackcrossing->Fill(trkcrossing);
    h_trackchi2ndf->Fill(chi2ndf);
    h_dcaxyorigin_phi->Fill(phi, dcapair_origin.first.first);
    h_dcaxyvtx_phi->Fill(phi, dcapair_vtx.first.first);
    h_dcazorigin_phi->Fill(phi, dcapair_origin.second.first);
    h_dcazvtx_phi->Fill(phi, dcapair_vtx.second.first);
    h_trackpt_inclusive->Fill(pt);
    if (charge > 0)
    {
      h_trackpt_pos->Fill(pt);
    }
    else
    {
      h_trackpt_neg->Fill(pt);
    }
  }

  h_ntrack_isfromvtx->SetBinContent(1, h_ntrack_isfromvtx->GetBinContent(1) + ntrack_isfromvtx.first);
  h_ntrack_isfromvtx->SetBinContent(2, h_ntrack_isfromvtx->GetBinContent(2) + ntrack_isfromvtx.second);

  h_ntrack_IsPosCharge->SetBinContent(1, h_ntrack_IsPosCharge->GetBinContent(1) + ntrack_isposcharge.first);
  h_ntrack_IsPosCharge->SetBinContent(2, h_ntrack_IsPosCharge->GetBinContent(2) + ntrack_isposcharge.second);

  // vertex
  h_nvertex->Fill(vertexmap->size());
  for (const auto &[key, vertex] : *vertexmap)
  {
    if (!vertex)
    {
      continue;
    }

    float vx = vertex->get_x();
    float vy = vertex->get_y();
    float vz = vertex->get_z();
    float vt = vertex->get_t0();
    float vchi2 = vertex->get_chisq();
    int vndof = vertex->get_ndof();
    int vcrossing = vertex->get_beam_crossing();

    // std::cout << "vertex (x,y,z,t,chi2,ndof,crossing)=(" << vx << "," << vy << "," << vz << "," << vt << "," << vchi2 << "," << vndof << "," << vcrossing << ")" << std::endl;

    h_vx->Fill(vx);
    h_vy->Fill(vy);
    h_vx_vy->Fill(vx, vy);
    h_vz->Fill(vz);
    h_vt->Fill(vt);
    h_vchi2dof->Fill(float(vchi2 / vndof));
    h_vcrossing->Fill(vcrossing);

    h_ntrackpervertex->Fill(vertex->size_tracks());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
std::vector<TrkrDefs::cluskey> SiliconSeedsQA::get_cluster_keys(SvtxTrack *track)
{
  std::vector<TrkrDefs::cluskey> out;
  for (const auto &seed : {track->get_silicon_seed(), track->get_tpc_seed()})
  {
    if (seed)
    {
      std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
    }
  }
  return out;
}

//____________________________________________________________________________..
int SiliconSeedsQA::EndRun(const int /*runnumber*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SiliconSeedsQA::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string SiliconSeedsQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void SiliconSeedsQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  {
    h_nmaps = new TH1F(std::string(getHistoPrefix() + "nmaps").c_str(), "MVTX clusters per track;Number of MVTX clusters per track;Entries", 7, -0.5, 6.5);
    hm->registerHisto(h_nmaps);
  }

  {
    h_nintt = new TH1F(std::string(getHistoPrefix() + "nintt").c_str(), "INTT clusters per track;Number of INTT clusters per track;Entries", 7, -0.5, 6.5);
    hm->registerHisto(h_nintt);
  }

  {
    h_nmaps_nintt = new TH2F(std::string(getHistoPrefix() + "nmaps_nintt").c_str(), "MVTX vs INTT clusters per track;Number of MVTX clusters per track;Number of INTT clusters per track;Entries", 7, -0.5, 6.5, 7, -0.5, 6.5);
    hm->registerHisto(h_nmaps_nintt);
  }

  {
    h_ntrack1d = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d").c_str(), "Number of reconstructed tracks;Number of silicon tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h_ntrack1d);
  }

  {
    h_ntrack = new TH2F(std::string(getHistoPrefix() + "nrecotracks").c_str(), "Number of reconstructed tracks;#eta;#phi [rad];Entries", 100, -1.1, 1.1, 300, -3.14159, 3.1459);
    hm->registerHisto(h_ntrack);
  }

  {
    h_avgnclus_eta_phi = new TProfile2D(std::string(getHistoPrefix() + "avgnclus_eta_phi").c_str(), "Average number of clusters per track;#eta;#phi [rad];Average number of clusters per track", 100, -1.1, 1.1, 300, -3.14159, 3.1459, 0, 10);
    hm->registerHisto(h_avgnclus_eta_phi);
  }

  {
    h_trackcrossing = new TH1F(std::string(getHistoPrefix() + "trackcrossing").c_str(), "Track beam bunch crossing;Track crossing;Entries", 110, -10, 100);
    hm->registerHisto(h_trackcrossing);
  }

  {
    h_trackchi2ndf = new TH1F(std::string(getHistoPrefix() + "trackchi2ndf").c_str(), "Track chi2/ndof;Track #chi2/ndof;Entries", 100, 0, 20);
    hm->registerHisto(h_trackchi2ndf);
  }

  {
    h_dcaxyorigin_phi = new TH2F(std::string(getHistoPrefix() + "dcaxyorigin_phi").c_str(), "DCA xy origin vs phi;#phi [rad];DCA_{xy} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 90, -3, 3);
    hm->registerHisto(h_dcaxyorigin_phi);
  }

  {
    h_dcaxyvtx_phi = new TH2F(std::string(getHistoPrefix() + "dcaxyvtx_phi").c_str(), "DCA xy vertex vs phi;#phi [rad];DCA_{xy} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 90, -3, 3);
    hm->registerHisto(h_dcaxyvtx_phi);
  }

  {
    h_dcazorigin_phi = new TH2F(std::string(getHistoPrefix() + "dcazorigin_phi").c_str(), "DCA z origin vs phi;#phi [rad];DCA_{z} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 160, -20, 20);
    hm->registerHisto(h_dcazorigin_phi);
  }

  {
    h_dcazvtx_phi = new TH2F(std::string(getHistoPrefix() + "dcazvtx_phi").c_str(), "DCA z vertex vs phi;#phi [rad];DCA_{z} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 160, -20, 20);
    hm->registerHisto(h_dcazvtx_phi);
  }

  {
    h_ntrack_isfromvtx = new TH1F(std::string(getHistoPrefix() + "ntrack_isfromvtx").c_str(), "Num of tracks associated to a vertex;Is track associated to a vertex;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h_ntrack_isfromvtx);
  }

  {
    h_trackpt_inclusive = new TH1F(std::string(getHistoPrefix() + "trackpt").c_str(), "Track p_{T};p_{T} [GeV/c];Entries", 100, 0, 10);
    hm->registerHisto(h_trackpt_inclusive);
  }

  {
    h_trackpt_pos = new TH1F(std::string(getHistoPrefix() + "trackpt_pos").c_str(), "Track p_{T} positive;Positive track p_{T} [GeV/c];Entries", 100, 0, 10);
    hm->registerHisto(h_trackpt_pos);
  }

  {
    h_trackpt_neg = new TH1F(std::string(getHistoPrefix() + "trackpt_neg").c_str(), "Track p_{T} negative;Negative track p_{T} [GeV/c];Entries", 100, 0, 10);
    hm->registerHisto(h_trackpt_neg);
  }

  {
    h_ntrack_IsPosCharge = new TH1F(std::string(getHistoPrefix() + "ntrack_IsPosCharge").c_str(), "Num of tracks with positive charge;Is track positive charged;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h_ntrack_IsPosCharge);
  }

  // vertex
  {
    h_nvertex = new TH1F(std::string(getHistoPrefix() + "nrecovertices").c_str(), "Num of reco vertices per event;Number of vertices;Entries", 20, 0, 20);
    hm->registerHisto(h_nvertex);
  }

  {
    h_vx = new TH1F(std::string(getHistoPrefix() + "vx").c_str(), "Vertex x;Vertex x [cm];Entries", 100, -2.5, 2.5);
    hm->registerHisto(h_vx);
  }

  {
    h_vy = new TH1F(std::string(getHistoPrefix() + "vy").c_str(), "Vertex y;Vertex y [cm];Entries", 100, -2.5, 2.5);
    hm->registerHisto(h_vy);
  }

  {
    h_vx_vy = new TH2F(std::string(getHistoPrefix() + "vx_vy").c_str(), "Vertex x vs y;Vertex x [cm];Vertex y [cm];Entries", 100, -2.5, 2.5, 100, -2.5, 2.5);
    hm->registerHisto(h_vx_vy);
  }

  {
    h_vz = new TH1F(std::string(getHistoPrefix() + "vz").c_str(), "Vertex z;Vertex z [cm];Entries", 50, -25, 25);
    hm->registerHisto(h_vz);
  }

  {
    h_vt = new TH1F(std::string(getHistoPrefix() + "vt").c_str(), "Vertex t;Vertex t [ns];Entries", 100, -1000, 20000);
    hm->registerHisto(h_vt);
  }

  {
    h_vcrossing = new TH1F(std::string(getHistoPrefix() + "vertexcrossing").c_str(), "Vertex beam bunch crossing;Vertex crossing;Entries", 100, -100, 300);
    hm->registerHisto(h_vcrossing);
  }

  {
    h_vchi2dof = new TH1F(std::string(getHistoPrefix() + "vertexchi2dof").c_str(), "Vertex chi2/ndof;Vertex #chi2/ndof;Entries", 100, 0, 20);
    hm->registerHisto(h_vchi2dof);
  }

  {
    h_ntrackpervertex = new TH1F(std::string(getHistoPrefix() + "ntrackspervertex").c_str(), "Num of tracks per vertex;Number of tracks per vertex;Entries", 50, 0, 50);
    hm->registerHisto(h_ntrackpervertex);
  }
}
