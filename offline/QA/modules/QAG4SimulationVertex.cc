#include "QAG4SimulationVertex.h"
#include "QAHistManagerDef.h"

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>
#include <g4eval/SvtxTruthEval.h>
#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxVertexEval.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h> // for cluskey

#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>
#include <TAxis.h>
#include <TDatabasePDG.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TParticlePDG.h>  // for TParticlePDG
#include <TString.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

using namespace std;

QAG4SimulationVertex::QAG4SimulationVertex(const std::string &name)
  : SubsysReco(name)
{
  cout << "QAG4SimulationVertex::QAG4SimulationVertex(const std::string &name) Calling ctor" << endl;
}

int QAG4SimulationVertex::InitRun(PHCompositeNode *topNode)
{
  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack.reset(new SvtxEvalStack(topNode));
    m_svtxEvalStack->set_strict(false);
    m_svtxEvalStack->set_verbosity(Verbosity() + 1);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationVertex::Init(PHCompositeNode *topNode)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h(nullptr);

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
               "ntracks Distribution;Number of Tracks;Count", 100, 0.5, 100.5);
  hm->registerHisto(h);

  // gntracks distribution histogram
  h = new TH1F(TString(get_histo_prefix()) + "gntracks",
               "gntracks Distribution;Number of gTracks;Count", 50, 0.5, 100.5);
  hm->registerHisto(h);

  // gntracksmaps distibution histogram
  h = new TH1F(TString(get_histo_prefix()) + "gntracksmaps",
               "gntracksmaps Distribution;Number of gTracksMaps;Count", 50, 0.5, 100.5);
  hm->registerHisto(h);

  // gvz distribution histogram
  h = new TH1F(TString(get_histo_prefix()) + "gvz",
               "gvz Distribution;gvz Position [cm]", 100, -10., 10.);
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
    cout << "QAG4SimulationVertex::process_event() entered" << endl;

  // load relevant nodes from NodeTree
  load_nodes(topNode);

  // histogram manager
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  
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

  SvtxVertexMap *vertexmap = nullptr;

  int n_recoSvtxVertex = 0;

  vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMapRefit");
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (vertexmap && truthinfo)
  {
    const auto prange = truthinfo->GetPrimaryParticleRange();
    map<int, unsigned int> embedvtxid_particle_count;
    map<int, unsigned int> embedvtxid_maps_particle_count;
    map<int, unsigned int> vertex_particle_count;
    for (auto iter = prange.first; iter != prange.second; ++iter)  // process all primary paricle
    {
      const int point_id = iter->second->get_vtx_id();
      int gembed = truthinfo->isEmbededVtx(iter->second->get_vtx_id());
      ++vertex_particle_count[point_id];
      ++embedvtxid_particle_count[gembed];
      PHG4Particle* g4particle = iter->second;
      
      if (false && gembed <= 0) continue;
      
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
        //cout<<__LINE__<<": " << _ievent <<": " <<gtrackID << ": " << layer <<": " <<g4cluster->get_id() <<endl;
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

    auto vrange = truthinfo->GetPrimaryVtxRange();
    map<int, bool> embedvtxid_found;
    map<int, int> embedvtxid_vertex_id;
    map<int, PHG4VtxPoint*> embedvtxid_vertex;
    for (auto iter = vrange.first; iter != vrange.second; ++iter)  // process all primary vertexes
    {
      const int point_id = iter->first;
      int gembed = truthinfo->isEmbededVtx(point_id);
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
    
    for (SvtxVertexMap::Iter iter = vertexmap->begin();
    iter != vertexmap->end();
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
           float gntracks = truthinfo->GetNumPrimaryVertexParticles();
           float gntracksmaps = NAN;

           h_ntracks->Fill(ntracks);

           if (point)
           {
             const int point_id = point->get_id();
             gvx = point->get_x();
             gvy = point->get_y();
             gvz = point->get_z();
             gvt = point->get_t();

             h_gvz->Fill(gvz);

             gembed = truthinfo->isEmbededVtx(point_id);
             gntracks = embedvtxid_particle_count[(int) gembed];

             h_gntracks->Fill(gntracks);

             if (embedvtxid_maps_particle_count[(int) gembed] > 0 && fabs(gvt) < 2000. && fabs(gvz) < 13.0)
               gntracksmaps = embedvtxid_maps_particle_count[(int) gembed];

             h_gntracksmaps->Fill(gntracksmaps);
           }
           float vx_res = vx - gvx;
           float vy_res = vy - gvy;
           float vz_res = vz - gvz;

           h_vxRes_gvz->Fill(gvz, vx_res);
           h_vyRes_gvz->Fill(gvz, vy_res);
           h_vzRes_gvz->Fill(gvz, vz_res);
         }
         h_recoSvtxVertex->Fill(n_recoSvtxVertex);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationVertex::load_nodes(PHCompositeNode *topNode)
{
  m_truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthContainer)
  {
    cout << "QAG4SimulationTracking::load_nodes - Fatal Error - "
         << "unable to find DST node "
         << "G4TruthInfo" << endl;
    assert(m_truthContainer);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

string
QAG4SimulationVertex::get_histo_prefix()
{
  return string("h_") + Name() + string("_");
}

