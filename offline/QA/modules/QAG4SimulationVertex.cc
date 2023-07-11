#include "QAG4SimulationVertex.h"

#include "QAHistManagerDef.h"

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxVertexEval.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/getClass.h>

#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>
#include <TVector3.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>  // for pair

QAG4SimulationVertex::QAG4SimulationVertex(const std::string &name)
  : SubsysReco(name)
{
}

int QAG4SimulationVertex::InitRun(PHCompositeNode *topNode)
{
  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack.reset(new SvtxEvalStack(topNode));
    m_svtxEvalStack->set_strict(false);
    m_svtxEvalStack->set_verbosity(Verbosity());
  }
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing " << m_trackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapName);
  if (!m_trackMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing " << m_vertexMapName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_trackMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing G4TruthInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationVertex::Init(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h(nullptr);

  h = new TH1D(TString(get_histo_prefix()) + "Normalization",  //
               "Normalization;Items;Count", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event");
  h->GetXaxis()->SetBinLabel(i++, "Vertex");
  h->GetXaxis()->SetBinLabel(i++, "MVTXTrackOnVertex");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  // Vertex resolution histograms as a funciton of gvz
  h = new TH2F(TString(get_histo_prefix()) + "vxRes_gvz",
               "Vertex Resolution at gvz (x);gvz [cm];<vx> - <gvx> [cm]", 100, -10., 10., 500, -0.005, 0.005);
  // QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "vyRes_gvz",
               "Vertex Resolution at gvz (y);gvz [cm];<vy> - <gvy> [cm]", 100, -10., 10., 500, -0.005, 0.005);
  // QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "vzRes_gvz",
               "Vertex Resolution at gvz (z);gvz [cm];<vz> - <gvz> [cm]", 100, -10., 10., 500, -0.005, 0.005);
  // QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  // ntracks distribution histogram
  h = new TH1F(TString(get_histo_prefix()) + "ntracks",
               "ntracks Distribution;Number of Tracks;Count", 200, 0.5, 200.5);
  hm->registerHisto(h);

  // ntracks distribution histogram with mvtx cuts
  h = new TH1F(TString(get_histo_prefix()) + "ntracks_cuts",
               "ntracks Distribution (#geq 2 MVTX);Number of Tracks;Count", 200, 0.5, 200.5);
  hm->registerHisto(h);

  // gntracks distribution histogram
  h = new TH1F(TString(get_histo_prefix()) + "gntracks",
               "gntracks Distribution;Number of gTracks;Count", 200, 0.5, 200.5);
  hm->registerHisto(h);

  // gntracksmaps distibution histogram
  h = new TH1F(TString(get_histo_prefix()) + "gntracksmaps",
               "gntracksmaps Distribution;Number of gTracksMaps;Count", 200, 0.5, 200.5);
  hm->registerHisto(h);

  // gvz distribution histogram
  h = new TH1F(TString(get_histo_prefix()) + "gvz",
               "gvz Distribution;gvz Position [cm]", 300, -15., 15.);
  hm->registerHisto(h);

  // Reco SvtxVertex count per event histogram
  h = new TH1F(TString(get_histo_prefix()) + "recoSvtxVertex",
               "SvtxVertex Count per Event;Number of SvtxVertex", 50, -0.5, 49.5);
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
}

void QAG4SimulationVertex::addEmbeddingID(int embeddingID)
{
  m_embeddingIDs.insert(embeddingID);
}

int QAG4SimulationVertex::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
    std::cout << "QAG4SimulationVertex::process_event() entered" << std::endl;

  // load relevant nodes from NodeTree
  load_nodes(topNode);

  // histogram manager
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto(
      get_histo_prefix() + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);
  ;

  if (m_svtxEvalStack)
    m_svtxEvalStack->next_event(topNode);
  /*
  SvtxTrackEval *trackeval = m_svtxEvalStack->get_track_eval();
  assert(trackeval);
  SvtxTruthEval *trutheval = m_svtxEvalStack->get_truth_eval();
  assert(trutheval);
  */
  SvtxVertexEval *vertexeval = m_svtxEvalStack->get_vertex_eval();
  SvtxClusterEval *clustereval = m_svtxEvalStack->get_cluster_eval();

  // Vertex resolution histograms
  TH2 *h_vxRes_gvz = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "vxRes_gvz"));
  assert(h_vxRes_gvz);
  TH2 *h_vyRes_gvz = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "vyRes_gvz"));
  assert(h_vyRes_gvz);
  TH2 *h_vzRes_gvz = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "vzRes_gvz"));
  assert(h_vzRes_gvz);

  // ntracks distribution histogram
  TH1 *h_ntracks = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "ntracks"));
  assert(h_ntracks);

  // ntracks distribution histogram with mvtx cuts
  TH1 *h_ntracks_cuts = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "ntracks_cuts"));
  assert(h_ntracks_cuts);

  // gntracks histogram
  TH1 *h_gntracks = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "gntracks"));
  assert(h_gntracks);

  // gntracksmaps histogram
  TH1 *h_gntracksmaps = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "gntracksmaps"));
  assert(h_gntracksmaps);

  // gvz histogram
  TH1 *h_gvz = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "gvz"));
  assert(h_gvz);

  // Reco SvtxVertex Histogram
  TH1 *h_recoSvtxVertex = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "recoSvtxVertex"));
  assert(h_recoSvtxVertex);

  int n_recoSvtxVertex = 0;
  if (m_vertexMap && m_truthInfo)
  {
    const auto prange = m_truthInfo->GetPrimaryParticleRange();
    std::map<int, unsigned int> embedvtxid_particle_count;
    std::map<int, unsigned int> embedvtxid_maps_particle_count;
    std::map<int, unsigned int> vertex_particle_count;
    for (auto iter = prange.first; iter != prange.second; ++iter)  // process all primary paricle
    {
      const int point_id = iter->second->get_vtx_id();
      int gembed = m_truthInfo->isEmbededVtx(iter->second->get_vtx_id());
      ++vertex_particle_count[point_id];
      ++embedvtxid_particle_count[gembed];
      PHG4Particle *g4particle = iter->second;

      if (m_checkembed && gembed <= m_embed_id_cut) continue;

      std::set<TrkrDefs::cluskey> g4clusters = clustereval->all_clusters_from(g4particle);

      unsigned int nglmaps = 0;

      int lmaps[_nlayers_maps + 1];

      if (_nlayers_maps > 0)
      {
        for (unsigned int i = 0; i < _nlayers_maps; i++)
        {
          lmaps[i] = 0;
        }
      }

      for (const TrkrDefs::cluskey g4cluster : g4clusters)
      {
        unsigned int layer = TrkrDefs::getLayer(g4cluster);
        // std::cout<<__LINE__<<": " << _ievent <<": " <<gtrackID << ": " << layer <<": " <<g4cluster->get_id() <<std::endl;
        if (_nlayers_maps > 0 && layer < _nlayers_maps)
        {
          lmaps[layer] = 1;
        }
      }
      if (_nlayers_maps > 0)
      {
        for (unsigned int i = 0; i < _nlayers_maps; i++)
        {
          nglmaps += lmaps[i];
        }
      }

      float gpx = g4particle->get_px();
      float gpy = g4particle->get_py();
      float gpz = g4particle->get_pz();
      float gpt = NAN;
      float geta = NAN;

      if (gpx != 0 && gpy != 0)
      {
        TVector3 gv(gpx, gpy, gpz);
        gpt = gv.Pt();
        geta = gv.Eta();
        //          gphi = gv.Phi();
      }
      if (nglmaps == 3 && fabs(geta) < 1.0 && gpt > 0.5)
        ++embedvtxid_maps_particle_count[gembed];
    }

    auto vrange = m_truthInfo->GetPrimaryVtxRange();
    std::map<int, bool> embedvtxid_found;
    std::map<int, int> embedvtxid_vertex_id;
    std::map<int, PHG4VtxPoint *> embedvtxid_vertex;
    for (auto iter = vrange.first; iter != vrange.second; ++iter)  // process all primary vertexes
    {
      const int point_id = iter->first;
      int gembed = m_truthInfo->isEmbededVtx(point_id);
      if (gembed <= 0) continue;

      auto search = embedvtxid_found.find(gembed);
      if (search != embedvtxid_found.end())
      {
        embedvtxid_vertex_id[gembed] = point_id;
        embedvtxid_vertex[gembed] = iter->second;
      }
      else
      {
        if (vertex_particle_count[embedvtxid_vertex_id[gembed]] < vertex_particle_count[point_id])
        {
          embedvtxid_vertex_id[gembed] = point_id;
          embedvtxid_vertex[gembed] = iter->second;
        }
      }
      embedvtxid_found[gembed] = false;
    }

    unsigned int ngembed = 0;
    for (std::map<int, bool>::iterator iter = embedvtxid_found.begin();
         iter != embedvtxid_found.end();
         ++iter)
    {
      if (iter->first >= 0 || iter->first != iter->first) continue;
      ++ngembed;
    }

    for (SvtxVertexMap::Iter iter = m_vertexMap->begin();
         iter != m_vertexMap->end();
         ++iter)
    {
      SvtxVertex *vertex = iter->second;
      ++n_recoSvtxVertex;

      PHG4VtxPoint *point = vertexeval->max_truth_point_by_ntracks(vertex);
      float vx = vertex->get_x();
      float vy = vertex->get_y();
      float vz = vertex->get_z();
      float ntracks = vertex->size_tracks();
      float gvx = NAN;
      float gvy = NAN;
      float gvz = NAN;
      float gvt = NAN;
      float gembed = NAN;
      float gntracks = m_truthInfo->GetNumPrimaryVertexParticles();
      float gntracksmaps = NAN;

      h_ntracks->Fill(ntracks);

      int ntracks_with_cuts = 0;
      for (SvtxVertex::TrackIter iter2 = vertex->begin_tracks();
           iter2 != vertex->end_tracks();
           ++iter2)
      {
        SvtxTrack *track = m_trackMap->get(*iter2);

        if (false)
        {
          assert(track);
        }
        else if (!track)
        {
          continue;
        }
        int MVTX_hits = 0;
        int INTT_hits = 0;
        int TPC_hits = 0;

        TrackSeed *siliconSeed = track->get_tpc_seed();
        TrackSeed *tpcSeed = track->get_silicon_seed();
        if (siliconSeed)
        {
          for (auto cluster_iter = siliconSeed->begin_cluster_keys();
               cluster_iter != siliconSeed->end_cluster_keys(); ++cluster_iter)
          {
            const auto &cluster_key = *cluster_iter;
            const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

            if (trackerID == TrkrDefs::mvtxId)
              ++MVTX_hits;
            else if (trackerID == TrkrDefs::inttId)
              ++INTT_hits;
            else
            {
              if (Verbosity())
                std::cout << "QAG4SimulationTracking::process_event - unkown tracker ID = " << trackerID << " from cluster " << cluster_key << std::endl;
            }
          }
        }
        for (auto cluster_iter = tpcSeed->begin_cluster_keys();
             cluster_iter != tpcSeed->end_cluster_keys(); ++cluster_iter)
        {
          const auto &cluster_key = *cluster_iter;
          const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

          if (trackerID == TrkrDefs::tpcId)
            ++TPC_hits;
        }
        if (MVTX_hits >= 2)
        {
          ++ntracks_with_cuts;
        }
      }
      h_ntracks_cuts->Fill(ntracks_with_cuts);

      if (point)
      {
        const int point_id = point->get_id();
        gvx = point->get_x();
        gvy = point->get_y();
        gvz = point->get_z();
        gvt = point->get_t();

        h_gvz->Fill(gvz);

        gembed = m_truthInfo->isEmbededVtx(point_id);
        gntracks = embedvtxid_particle_count[(int) gembed];

        h_gntracks->Fill(gntracks);

        if (embedvtxid_maps_particle_count[(int) gembed] > 0 && fabs(gvt) < 2000. && fabs(gvz) < 13.0)
          gntracksmaps = embedvtxid_maps_particle_count[(int) gembed];

        h_gntracksmaps->Fill(gntracksmaps);
        h_norm->Fill("MVTXTrackOnVertex", gntracksmaps);
      }
      float vx_res = vx - gvx;
      float vy_res = vy - gvy;
      float vz_res = vz - gvz;

      h_vxRes_gvz->Fill(gvz, vx_res);
      h_vyRes_gvz->Fill(gvz, vy_res);
      h_vzRes_gvz->Fill(gvz, vz_res);
    }
    h_recoSvtxVertex->Fill(n_recoSvtxVertex);
    h_norm->Fill("Vertex", n_recoSvtxVertex);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationVertex::load_nodes(PHCompositeNode *topNode)
{
  m_truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthContainer)
  {
    std::cout << "QAG4SimulationTracking::load_nodes - Fatal Error - "
              << "unable to find DST node "
              << "G4TruthInfo" << std::endl;
    assert(m_truthContainer);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string
QAG4SimulationVertex::get_histo_prefix()
{
  return std::string("h_") + Name() + std::string("_") + m_vertexMapName + std::string("_");
}
