#include "TpcSeedsQA.h"

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>

//____________________________________________________________________________..
TpcSeedsQA::TpcSeedsQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TpcSeedsQA::InitRun(PHCompositeNode *topNode)
{
  createHistos();

  clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);

  if (!trackmap or !clustermap or !geometry or !vertexmap)
  {
    std::cout << PHWHERE << "Missing node(s), can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // tracks with TPC clusters/tracklets
  h_ntrack1d = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d").c_str()));
  h_ntrack1d_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_pos").c_str()));
  h_ntrack1d_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_neg").c_str()));
  h_ntrack1d_ptg1 = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_ptg1").c_str()));
  h_ntrack1d_ptg1_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_ptg1_pos").c_str()));
  h_ntrack1d_ptg1_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d_ptg1_neg").c_str()));
  h_pt = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "pt").c_str()));
  h_pt_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "pt_pos").c_str()));
  h_pt_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "pt_neg").c_str()));
  h_ntrack_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks_pos").c_str()));
  h_ntrack_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks_neg").c_str()));

  h_ntpc_fullpt_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_fullpt_pos").c_str()));
  h_ntpc_fullpt_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_fullpt_neg").c_str()));
  h_ntpc_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_pos").c_str()));
  h_ntpc_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_neg").c_str()));
  h_ntpc_quality_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_quality_pos").c_str()));
  h_ntpc_quality_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpc_quality_neg").c_str()));
  h_ntpot_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpot_pos").c_str()));
  h_ntpot_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntpot_neg").c_str()));
  h_avgnclus_eta_phi_pos = dynamic_cast<TProfile2D *>(hm->getHisto(std::string(getHistoPrefix() + "avgnclus_eta_phi_pos").c_str()));
  h_avgnclus_eta_phi_neg = dynamic_cast<TProfile2D *>(hm->getHisto(std::string(getHistoPrefix() + "avgnclus_eta_phi_neg").c_str()));
  //h_trackcrossing_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "trackcrossing_pos").c_str()));
  //h_trackcrossing_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "trackcrossing_neg").c_str()));
  h_dcaxyorigin_phi_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyorigin_phi_pos").c_str()));
  h_dcaxyorigin_phi_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyorigin_phi_neg").c_str()));
  h_dcaxyvtx_phi_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyvtx_phi_pos").c_str()));
  h_dcaxyvtx_phi_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcaxyvtx_phi_neg").c_str()));
  h_dcazorigin_phi_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcazorigin_phi_pos").c_str()));
  h_dcazorigin_phi_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcazorigin_phi_neg").c_str()));
  h_dcazvtx_phi_pos = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcazvtx_phi_pos").c_str()));
  h_dcazvtx_phi_neg = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "dcazvtx_phi_neg").c_str()));
  h_ntrack_isfromvtx_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntrack_isfromvtx_pos").c_str()));
  h_ntrack_isfromvtx_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntrack_isfromvtx_neg").c_str()));
  h_cluster_phisize1_fraction_pos = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "cluster_phisize1_fraction_pos").c_str()));
  h_cluster_phisize1_fraction_neg = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "cluster_phisize1_fraction_neg").c_str()));

  // vertex
  h_nvertex = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecovertices").c_str()));
  h_vx = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vx").c_str()));
  h_vy = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vy").c_str()));
  h_vx_vy = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "vx_vy").c_str()));
  h_vz = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vz").c_str()));
  h_vt = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vt").c_str()));
  //h_vcrossing = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexcrossing").c_str()));
  h_vchi2dof = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexchi2dof").c_str()));
  h_ntrackpervertex = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntrackspervertex").c_str()));

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcSeedsQA::process_event(PHCompositeNode * /*unused*/)
{
  h_ntrack1d->Fill(trackmap->size());

  std::pair<int, int> ntrack_isfromvtx_pos;  // first: number of tracks not associated to a vertex, second: number of tracks associated to a vertex
  std::pair<int, int> ntrack_isfromvtx_neg;  // first: number of tracks not associated to a vertex, second: number of tracks associated to a vertex

  int ntrack1d_pos = 0;
  int ntrack1d_neg = 0;
  int ntrack1d_ptg1_pos = 0;
  int ntrack1d_ptg1_neg = 0;

  for (const auto &[key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }

    int charge = track->get_charge();
    float quality = track->get_quality();
    float pt = track->get_pt();

    h_pt->Fill(pt);
    if (charge == 1)
    {
      ntrack1d_pos++;
      if (pt>1) ntrack1d_ptg1_pos++;
      h_pt_pos->Fill(pt);
    }
    else if (charge == -1)
    {
      ntrack1d_neg++;
      if (pt>1) ntrack1d_ptg1_neg++;
      h_pt_neg->Fill(pt);
    }

    auto ckeys = get_cluster_keys(track);
    std::vector<Acts::Vector3> cluspos;
    TrackFitUtils::getTrackletClusters(geometry, clustermap, cluspos, ckeys);
    float eta = track->get_eta();
    float phi = track->get_phi();

    //int trkcrossing = track->get_crossing();

    int nmaps = 0;
    int nintt = 0;
    int ntpc = 0;
    int ntpc_phisize1 = 0;
    int nmms = 0;

    for (auto &ckey : ckeys)
    {
      if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::tpcId)
      {
        TrkrCluster *cluster = clustermap->findCluster(ckey);
        if (cluster->getPhiSize() == 1)
        {
          ntpc_phisize1++;
        }
      }
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

    Acts::Vector3 zero = Acts::Vector3::Zero();
    auto dcapair_origin = TrackAnalysisUtils::get_dca(track, zero);

    auto trackvtx = vertexmap->get(track->get_vertex_id());
    if (!trackvtx)
    {
      // std::cout << "No vertex found for track " << track->get_id() << std::std::endl;
      if (charge == 1)
      {
        ntrack_isfromvtx_pos.first++;
      }
      else if (charge == -1)
      {
        ntrack_isfromvtx_neg.first++;
      }
    }
    else
    {
      Acts::Vector3 track_vtx(trackvtx->get_x(), trackvtx->get_y(), trackvtx->get_z());
      auto dcapair_vtx = TrackAnalysisUtils::get_dca(track, track_vtx);
      if (charge == 1)
      {
        ntrack_isfromvtx_pos.second++;
        h_dcaxyvtx_phi_pos->Fill(phi, dcapair_vtx.first.first);
        h_dcazvtx_phi_pos->Fill(phi, dcapair_vtx.second.first);
      }
      else if (charge == -1)
      {
        ntrack_isfromvtx_neg.second++;
        h_dcaxyvtx_phi_neg->Fill(phi, dcapair_vtx.first.first);
        h_dcazvtx_phi_neg->Fill(phi, dcapair_vtx.second.first);
      }
    }

    if (charge == 1)
    {
      h_ntpc_fullpt_pos->Fill(ntpc);
      h_dcaxyorigin_phi_pos->Fill(phi, dcapair_origin.first.first);
      h_dcazorigin_phi_pos->Fill(phi, dcapair_origin.second.first);
      if (pt>1)
      {
        h_ntrack_pos->Fill(eta, phi);
        h_ntpc_pos->Fill(ntpc);
        h_ntpot_pos->Fill(nmms);
        h_ntpc_quality_pos->Fill(ntpc,quality);
        h_avgnclus_eta_phi_pos->Fill(eta, phi, ntpc);
        //h_trackcrossing_pos->Fill(trkcrossing);
        h_cluster_phisize1_fraction_pos->Fill((double) ntpc_phisize1 / (double) ntpc);
      }
    }
    else if (charge == -1)
    {
      h_ntpc_fullpt_neg->Fill(ntpc);
      h_dcaxyorigin_phi_neg->Fill(phi, dcapair_origin.first.first);
      h_dcazorigin_phi_neg->Fill(phi, dcapair_origin.second.first);
      if (pt>1)
      {
        h_ntrack_neg->Fill(eta, phi);
        h_ntpc_neg->Fill(ntpc);
        h_ntpot_neg->Fill(nmms);
        h_ntpc_quality_neg->Fill(ntpc,quality);
        h_avgnclus_eta_phi_neg->Fill(eta, phi, ntpc);
        //h_trackcrossing_neg->Fill(trkcrossing);
        h_cluster_phisize1_fraction_neg->Fill((double) ntpc_phisize1 / (double) ntpc);
      }
    }
  }
  h_ntrack1d_pos->Fill(ntrack1d_pos);
  h_ntrack1d_neg->Fill(ntrack1d_neg);
  h_ntrack1d_ptg1_pos->Fill(ntrack1d_ptg1_pos);
  h_ntrack1d_ptg1_neg->Fill(ntrack1d_ptg1_neg);
  h_ntrack1d_ptg1->Fill(ntrack1d_ptg1_pos + ntrack1d_ptg1_neg);

  h_ntrack_isfromvtx_pos->SetBinContent(1, h_ntrack_isfromvtx_pos->GetBinContent(1) + ntrack_isfromvtx_pos.first);
  h_ntrack_isfromvtx_pos->SetBinContent(2, h_ntrack_isfromvtx_pos->GetBinContent(2) + ntrack_isfromvtx_pos.second);
  h_ntrack_isfromvtx_neg->SetBinContent(1, h_ntrack_isfromvtx_neg->GetBinContent(1) + ntrack_isfromvtx_neg.first);
  h_ntrack_isfromvtx_neg->SetBinContent(2, h_ntrack_isfromvtx_neg->GetBinContent(2) + ntrack_isfromvtx_neg.second);

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
    //int vcrossing = vertex->get_beam_crossing();

    // std::cout << "vertex (x,y,z,t,chi2,ndof,crossing)=(" << vx << "," << vy << "," << vz << "," << vt << "," << vchi2 << "," << vndof << "," << vcrossing << ")" << std::endl;

    h_vx->Fill(vx);
    h_vy->Fill(vy);
    h_vx_vy->Fill(vx, vy);
    h_vz->Fill(vz);
    h_vt->Fill(vt);
    h_vchi2dof->Fill(float(vchi2 / vndof));
    //h_vcrossing->Fill(vcrossing);

    h_ntrackpervertex->Fill(vertex->size_tracks());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::vector<TrkrDefs::cluskey> TpcSeedsQA::get_cluster_keys(SvtxTrack *track)
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
int TpcSeedsQA::EndRun(const int /*runnumber*/)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcSeedsQA::End(PHCompositeNode * /*unused*/) { return Fun4AllReturnCodes::EVENT_OK; }

std::string TpcSeedsQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

void TpcSeedsQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpc_fullpt_pos").c_str(), "TPC clusters per positive track;Number of TPC clusters per positive track;Entries", 55, -0.5, 54.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpc_fullpt_neg").c_str(), "TPC clusters per negative track;Number of TPC clusters per negative track;Entries", 55, -0.5, 54.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpc_pos").c_str(), "TPC clusters per positive track (pT>1GeV);Number of TPC clusters per positive track;Entries", 55, -0.5, 54.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpc_neg").c_str(), "TPC clusters per negative track (pT>1GeV);Number of TPC clusters per negative track;Entries", 55, -0.5, 54.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpot_pos").c_str(), "TPOT clusters per positive track (pT>1GeV);Number of TPOT clusters per positive track;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntpot_neg").c_str(), "TPOT clusters per negative track (pT>1GeV);Number of TPOT clusters per negative track;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "ntpc_quality_pos").c_str(), "Number of TPC clusters per positive track (pT>1GeV);Number of TPC clusters per positive track;Quality", 55, -0.5, 54.5, 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "ntpc_quality_neg").c_str(), "Number of TPC clusters per negative track (pT>1GeV);Number of TPC clusters per negative track;Quality", 55, -0.5, 54.5, 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d").c_str(), "Number of reconstructed tracks;Number of TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_pos").c_str(), "Number of reconstructed positive tracks;Number of positive TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_neg").c_str(), "Number of reconstructed negative tracks;Number of negative TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_ptg1").c_str(), "Number of reconstructed tracks (pT>1GeV);Number of TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_ptg1_pos").c_str(), "Number of reconstructed positive tracks (pT>1GeV);Number of positive TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d_ptg1_neg").c_str(), "Number of reconstructed negative tracks (pT>1GeV);Number of negative TPC tracklets;Entries", 50, 0, 200);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "pt").c_str(), "p_{T} distribution of reconstructed tracks;Track p_{T};Entries", 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "pt_pos").c_str(), "p_{T} distribution of reconstructed positive tracks;Track p_{T};Entries", 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "pt_neg").c_str(), "p_{T} distribution of reconstructed negative tracks;Track p_{T};Entries", 100, 0, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "nrecotracks_pos").c_str(), "Number of reconstructed positive tracks (pT>1GeV);#eta;#phi [rad];Entries", 100, -1.1, 1.1, 300, -3.14159, 3.1459);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "nrecotracks_neg").c_str(), "Number of reconstructed negative tracks (pT>1GeV);#eta;#phi [rad];Entries", 100, -1.1, 1.1, 300, -3.14159, 3.1459);
    hm->registerHisto(h);
  }

  {
    auto h = new TProfile2D(std::string(getHistoPrefix() + "avgnclus_eta_phi_pos").c_str(), "Average number of clusters per positive track (pT>1GeV);#eta;#phi [rad];Average number of clusters per positive track", 100, -1.1, 1.1, 300, -3.14159, 3.1459, 0, 55);
    hm->registerHisto(h);
  }

  {
    auto h = new TProfile2D(std::string(getHistoPrefix() + "avgnclus_eta_phi_neg").c_str(), "Average number of clusters per negative track (pT>1GeV);#eta;#phi [rad];Average number of clusters per negative track", 100, -1.1, 1.1, 300, -3.14159, 3.1459, 0, 55);
    hm->registerHisto(h);
  }

//  {
//    auto h = new TH1F(std::string(getHistoPrefix() + "trackcrossing_pos").c_str(), "Positive track beam bunch crossing (pT>1GeV);Positive track crossing;Entries", 100, -100, 300);
//    hm->registerHisto(h);
//  }

//  {
//    auto h = new TH1F(std::string(getHistoPrefix() + "trackcrossing_neg").c_str(), "Negative track beam bunch crossing (pT>1GeV);Negative track crossing;Entries", 100, -100, 300);
//    hm->registerHisto(h);
//  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyorigin_phi_pos").c_str(), "DCA xy origin vs phi for positive track;#phi [rad];DCA_{xy} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 90, -3, 3);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyorigin_phi_neg").c_str(), "DCA xy origin vs phi for negative track;#phi [rad];DCA_{xy} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 90, -3, 3);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyvtx_phi_pos").c_str(), "DCA xy vertex vs phi for positive track;#phi [rad];DCA_{xy} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 90, -3, 3);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcaxyvtx_phi_neg").c_str(), "DCA xy vertex vs phi for negative track;#phi [rad];DCA_{xy} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 90, -3, 3);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcazorigin_phi_pos").c_str(), "DCA z origin vs phi for positive track;#phi [rad];DCA_{z} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcazorigin_phi_neg").c_str(), "DCA z origin vs phi for negative track;#phi [rad];DCA_{z} wrt origin [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcazvtx_phi_pos").c_str(), "DCA z vertex vs phi for positive track;#phi [rad];DCA_{z} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "dcazvtx_phi_neg").c_str(), "DCA z vertex vs phi for negative track;#phi [rad];DCA_{z} wrt vertex [cm];Entries", 300, -3.14159, 3.1459, 100, -10, 10);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntrack_isfromvtx_pos").c_str(), "Num of positive tracks associated to a vertex;Is track associated to a vertex;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntrack_isfromvtx_neg").c_str(), "Num of negative tracks associated to a vertex;Is track associated to a vertex;Entries", 2, -0.5, 1.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "cluster_phisize1_fraction_pos").c_str(), "Fraction of TPC clusters per positive track with phi size of 1 (pT>1GeV);Fraction of TPC clusters phi size of 1;Entries", 100, 0, 1);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "cluster_phisize1_fraction_neg").c_str(), "Fraction of TPC clusters per negative track with phi size of 1 (pT>1GeV);Fraction of TPC clusters phi size of 1;Entries", 100, 0, 1);
    hm->registerHisto(h);
  }

  // vertex
  {
    auto h = new TH1F(std::string(getHistoPrefix() + "nrecovertices").c_str(), "Num of reco vertices per event;Number of vertices;Entries", 20, 0, 20);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vx").c_str(), "Vertex x;Vertex x [cm];Entries", 100, -2.5, 2.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vy").c_str(), "Vertex y;Vertex y [cm];Entries", 100, -2.5, 2.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(std::string(getHistoPrefix() + "vx_vy").c_str(), "Vertex x vs y;Vertex x [cm];Vertex y [cm];Entries", 100, -2.5, 2.5, 100, -2.5, 2.5);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vz").c_str(), "Vertex z;Vertex z [cm];Entries", 50, -25, 25);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vt").c_str(), "Vertex t;Vertex t [ns];Entries", 100, -1000, 20000);
    hm->registerHisto(h);
  }

//  {
//    auto h = new TH1F(std::string(getHistoPrefix() + "vertexcrossing").c_str(), "Vertex beam bunch crossing;Vertex crossing;Entries", 100, -100, 300);
//    hm->registerHisto(h);
//  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "vertexchi2dof").c_str(), "Vertex chi2/ndof;Vertex #chi2/ndof;Entries", 100, 0, 20);
    hm->registerHisto(h);
  }

  {
    auto h = new TH1F(std::string(getHistoPrefix() + "ntrackspervertex").c_str(), "Num of tracks per vertex;Number of tracks per vertex;Entries", 50, 0, 50);
    hm->registerHisto(h);
  }
}
