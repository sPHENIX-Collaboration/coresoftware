
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
#include <trackbase/TrackFitUtils.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

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
    auto clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
    auto geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    auto trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
    auto vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);

    if (!trackmap or !clustermap or !geometry or !vertexmap)
    {
        std::cout << PHWHERE << "Missing node(s), can't continue" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    auto hm = QAHistManagerDef::getHistoManager();
    assert(hm);

    // tracklets (tracks with MVTX+INTT clusters)
    auto h_ntrack1d = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks1d").c_str()));
    auto h_ntrack = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecotracks").c_str()));
    auto h_nmaps = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nmaps").c_str()));
    auto h_nintt = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nintt").c_str()));
    auto h_nmaps_nintt = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nmaps_nintt").c_str()));
    auto h_avgnclus_phi_eta = dynamic_cast<TProfile2D *>(hm->getHisto(std::string(getHistoPrefix() + "avgnclus_phi_eta").c_str()));

    // vertex
    auto h_nvertex = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecovertices").c_str()));
    auto h_vx = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vx").c_str()));
    auto h_vy = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vy").c_str()));
    auto h_vx_vy = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "vx_vy").c_str()));
    auto h_vz = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vz").c_str()));
    auto h_vt = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vt").c_str()));
    auto h_vcrossing = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexcrossing").c_str()));
    auto h_vchi2 = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexchi2").c_str()));
    auto h_vndof = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexndof").c_str()));
    auto h_ntrackpervertex = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntrackspervertex").c_str()));

    std::cout << "Number of tracks: " << trackmap->size() << std::endl;
    h_ntrack1d->Fill(trackmap->size());

    for (const auto &[key, track] : *trackmap)
    {
        if (!track)
        {
            continue;
        }

        auto ckeys = get_cluster_keys(track);
        std::vector<Acts::Vector3> cluspos;
        TrackFitUtils::getTrackletClusters(geometry, clustermap, cluspos, ckeys);

        float px = track->get_px();
        float py = track->get_py();
        float pz = track->get_pz();
        float pt = std::sqrt(QAG4Util::square(px) + QAG4Util::square(py));
        float eta = std::atanh(pz / std::sqrt(QAG4Util::square(pt) + QAG4Util::square(pz)));
        float phi = std::atan2(py, px);

        int nmaps = 0;
        int nintt = 0;
        int ntpc = 0;
        int nmms = 0;

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
            case TrkrDefs::tpcId:
                ntpc++;
                break;
            case TrkrDefs::micromegasId:
                nmms++;
                break;
            }
        }

        h_ntrack->Fill(phi, eta);
        h_nmaps->Fill(nmaps);
        h_nintt->Fill(nintt);
        h_nmaps_nintt->Fill(nmaps, nintt);
        h_avgnclus_phi_eta->Fill(phi, eta, nmaps + nintt);
    }

    // vertex
    m_vertices += vertexmap->size();
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

        std::cout << "vertex (vx,vy,vz)=(" << vx << "," << vy << "," << vz << ")" << std::endl;

        h_vx->Fill(vx);
        h_vy->Fill(vy);
        h_vx_vy->Fill(vx, vy);
        h_vz->Fill(vz);
        h_vt->Fill(vt);
        h_vchi2->Fill(vchi2);
        h_vndof->Fill(vndof);
        h_vcrossing->Fill(vcrossing);

        h_ntrackpervertex->Fill(vertex->size_tracks());
    }

    m_event++;
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
    auto hm = QAHistManagerDef::getHistoManager();
    assert(hm);

    return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SiliconSeedsQA::End(PHCompositeNode * /*unused*/) { return Fun4AllReturnCodes::EVENT_OK; }

std::string SiliconSeedsQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

void SiliconSeedsQA::createHistos()
{
    auto hm = QAHistManagerDef::getHistoManager();
    assert(hm);

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "nmaps").c_str(), "MVTX clusters per track", 11, -0.5, 10.5);
        h->GetXaxis()->SetTitle("Number of MVTX clusters");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "nintt").c_str(), "INTT clusters per track", 11, -0.5, 10.5);
        h->GetXaxis()->SetTitle("Number of INTT clusters");
        hm->registerHisto(h);
    }

    {
        auto h = new TH2F(std::string(getHistoPrefix() + "nmaps_nintt").c_str(), "MVTX vs INTT clusters per track", 6, -0.5, 5.5, 6, -0.5, 5.5);
        h->GetXaxis()->SetTitle("Number of MVTX clusters");
        h->GetYaxis()->SetTitle("Number of INTT clusters");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "nrecotracks1d").c_str(), "Num reconstructed tracks", 50, 0, 200);
        h->GetXaxis()->SetTitle("Number of tracklets");
        hm->registerHisto(h);
    }

    {
        auto h = new TH2F(std::string(getHistoPrefix() + "nrecotracks").c_str(), "Num reconstructed tracks", 300, -3.14159, 3.1459, 100, -1.1, 1.1);
        h->GetXaxis()->SetTitle("#phi [rad]");
        h->GetYaxis()->SetTitle("#eta");
        hm->registerHisto(h);
    }

    {
        auto h = new TProfile2D(std::string(getHistoPrefix() + "avgnclus_phi_eta").c_str(), "Average number of clusters per track", 300, -3.14159, 3.1459, 100, -1.1, 1.1, 0, 10);
        h->GetXaxis()->SetTitle("#phi [rad]");
        h->GetYaxis()->SetTitle("#eta");
        hm->registerHisto(h);
    }

    // vertex
    {
        auto h = new TH1F(std::string(getHistoPrefix() + "nrecovertices").c_str(), "Num of reco vertices per event", 30, 0, 30);
        h->GetXaxis()->SetTitle("nVertices");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "vx").c_str(), "Vertex x", 100, -2.5, 2.5);
        h->GetXaxis()->SetTitle("vx [cm]");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "vy").c_str(), "Vertex y", 100, -2.5, 2.5);
        h->GetXaxis()->SetTitle("vy [cm]");
        hm->registerHisto(h);
    }

    {
        auto h = new TH2F(std::string(getHistoPrefix() + "vx_vy").c_str(), "Vertex x vs y", 100, -2.5, 2.5, 100, -2.5, 2.5);
        h->GetXaxis()->SetTitle("vx [cm]");
        h->GetYaxis()->SetTitle("vy [cm]");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "vz").c_str(), "Vertex z", 100, -25, 25);
        h->GetXaxis()->SetTitle("vz [cm]");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "vt").c_str(), "Vertex t", 100, -1000, 20000);
        h->GetXaxis()->SetTitle("vt [ns]");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "vertexcrossing").c_str(), "Vertex beam bunch crossing", 100, -100, 300);
        h->GetXaxis()->SetTitle("vertex crossing");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "vertexchi2").c_str(), "Vertex chi2", 100, 0, 10000);
        h->GetXaxis()->SetTitle("vertex #chi2");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "vertexndof").c_str(), "Vertex ndof", 50, 0, 50);
        h->GetXaxis()->SetTitle("vertex ndof");
        hm->registerHisto(h);
    }

    {
        auto h = new TH1F(std::string(getHistoPrefix() + "ntrackspervertex").c_str(), "Num of tracks per vertex", 50, 0, 200);
        h->GetXaxis()->SetTitle("nTracks per vertex");
        hm->registerHisto(h);
    }
}
