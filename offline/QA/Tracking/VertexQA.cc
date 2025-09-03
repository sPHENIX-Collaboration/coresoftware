
#include "VertexQA.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TH2.h>
#include <cassert>
//____________________________________________________________________________..
VertexQA::VertexQA(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int VertexQA::InitRun(PHCompositeNode * /*unused*/)
{
  createHistos();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int VertexQA::process_event(PHCompositeNode *topNode)
{
  // auto trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto *vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);
  if (!vertexmap)
  {
    std::cout << PHWHERE << "Missing node(s), can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  auto *h_nvertex = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "nrecovertices")));
  auto *h_vx = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vx")));
  auto *h_vy = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vy")));
  auto *h_vz = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vz")));
  auto *h_vt = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vt")));
  auto *h_vcrossing = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexcrossing")));
  auto *h_vchi2 = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexchi2")));
  auto *h_vndof = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "vertexndof")));
  auto *h_ntrackpervertex = dynamic_cast<TH1 *>(hm->getHisto(std::string(getHistoPrefix() + "ntrackspervertex")));

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

    h_vx->Fill(vx);
    h_vy->Fill(vy);
    h_vz->Fill(vz);
    h_vt->Fill(vt);
    h_vchi2->Fill(vchi2);
    h_vndof->Fill(vndof);
    h_vcrossing->Fill(vcrossing);

    h_ntrackpervertex->Fill(vertex->size_tracks());
    // loop over all tracks on vertex
    // for (Vertex::TrackIter iter = vertex->begin_tracks();
    //     iter != vertex->end_tracks();
    //     ++iter)
    //{
    //  SvtxTrack* track = trackmap->get(*iter);

    //  if (!track)
    //  {
    //    continue;
    //  }

    //  float px = track->get_px();
    //  float py = track->get_py();
    //  float pz = track->get_pz();
    //  float pt = std::sqrt(QAG4Util::square(px) + QAG4Util::square(py));
    //  float eta = std::atanh(pz / std::sqrt(QAG4Util::square(pt) + QAG4Util::square(pz)));
    //  float phi = std::atan2(py, px);
    //  int tcrossing = track->get_crossing();
    //}
  }

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int VertexQA::EndRun(const int runnumber)
{
  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH2 *h_verticesperevent = dynamic_cast<TH2 *>(hm->getHisto(std::string(getHistoPrefix() + "nverticesperrun")));
  // NOLINTNEXTLINE(bugprone-integer-division)
  h_verticesperevent->Fill(runnumber, m_vertices / m_event);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int VertexQA::End(PHCompositeNode * /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string VertexQA::getHistoPrefix() const
{
  return std::string("h_") + Name() + std::string("_");
}

void VertexQA::createHistos()
{
  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  {
    auto *h = new TH1F(std::string(getHistoPrefix() + "nrecovertices").c_str(),
                       "Num of reco vertices per event", 200, 0, 200);
    h->GetXaxis()->SetTitle("nVertices");
    hm->registerHisto(h);
  }
  {
    auto *h = new TH1F(std::string(getHistoPrefix() + "vx").c_str(),
                       "Vertex x", 100, -0.1, 0.1);
    h->GetXaxis()->SetTitle("vx [cm]");
    hm->registerHisto(h);
  }
  {
    auto *h = new TH1F(std::string(getHistoPrefix() + "vy").c_str(),
                       "Vertex y", 100, -0.1, 0.1);
    h->GetXaxis()->SetTitle("vy [cm]");
    hm->registerHisto(h);
  }
  {
    auto *h = new TH1F(std::string(getHistoPrefix() + "vz").c_str(),
                       "Vertex z", 100, -15, 15);
    h->GetXaxis()->SetTitle("vz [cm]");
    hm->registerHisto(h);
  }
  {
    auto *h = new TH1F(std::string(getHistoPrefix() + "vt").c_str(),
                       "Vertex t", 100, -1000, 20000);
    h->GetXaxis()->SetTitle("vt [ns]");
    hm->registerHisto(h);
  }
  {
    auto *h = new TH1F(std::string(getHistoPrefix() + "vertexcrossing").c_str(),
                       "Vertex beam bunch crossing", 100, -100, 300);
    h->GetXaxis()->SetTitle("vertex crossing");
    hm->registerHisto(h);
  }
  {
    auto *h = new TH1F(std::string(getHistoPrefix() + "vertexchi2").c_str(),
                       "Vertex chi2", 100, 0, 10000);
    h->GetXaxis()->SetTitle("vertex #chi2");
    hm->registerHisto(h);
  }
  {
    auto *h = new TH1F(std::string(getHistoPrefix() + "vertexndof").c_str(),
                       "Vertex ndof", 50, 0, 50);
    h->GetXaxis()->SetTitle("vertex ndof");
    hm->registerHisto(h);
  }
  {
    auto *h = new TH1F(std::string(getHistoPrefix() + "ntrackspervertex").c_str(),
                       "Num of tracks per vertex", 20, 0, 20);
    h->GetXaxis()->SetTitle("nTracks per vertex");
    hm->registerHisto(h);
  }

  {
    auto *h = new TH2F(std::string(getHistoPrefix() + "nverticesperrun").c_str(),
                       "Num reconstructed vertices per run", m_runbins, m_beginRun, m_endRun, 100, 0, 1);
    h->GetYaxis()->SetTitle("N_{vertives}/event");
    h->GetXaxis()->SetTitle("Run number");
    hm->registerHisto(h);
  }
}
