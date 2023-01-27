#include "QAG4SimulationTracking.h"
#include "QAHistManagerDef.h"

#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>  // for SvtxTrackEval
#include <g4eval/SvtxTruthEval.h>  // for SvtxTruthEval

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>  // for cluskey, getLayer
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TAxis.h>
#include <TDatabasePDG.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TParticlePDG.h>  // for TParticlePDG
#include <TString.h>
#include <TVector3.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>      // for map
#include <utility>  // for pair

QAG4SimulationTracking::QAG4SimulationTracking(const std::string &name)
  : SubsysReco(name)
{
}

int QAG4SimulationTracking::InitRun(PHCompositeNode *topNode)
{
  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack.reset(new SvtxEvalStack(topNode));
    m_svtxEvalStack->set_strict(false);
    m_svtxEvalStack->set_verbosity(Verbosity());
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  if (!m_trackMap or !m_vertexMap)
  {
    std::cout << PHWHERE << " missing track related container(s). Quitting"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationTracking::Init(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // reco pT / gen pT histogram
  TH1 *h(nullptr);

  h = new TH1F(TString(get_histo_prefix()) + "pTRecoGenRatio",
               ";Reco p_{T}/Truth p_{T}", 500, 0, 2);
  hm->registerHisto(h);

  h = new TH2F(TString(get_histo_prefix()) + "pTRecoGenRatio_pTGen",
               ";Truth p_{T} [GeV/c];Reco p_{T}/Truth p_{T}", 200, 0, 50, 500, 0, 2);
  //  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  // reco track w/ truth-track matched vs reco pT histograms
  h = new TH1F(TString(get_histo_prefix()) + "nGen_pTReco",
               "Gen tracks at reco p_{T}; Reco p_{T} [GeV/c]", 200, 0.1, 50.5);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "nReco_pTReco",
               ";Gen p_{T} [GeV/c];Track count / bin", 200, 0.1, 50.5);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "pTRecoTruthMatchedRatio_pTReco",
               ";Reco p_{T} [GeV/c];Matched p_{T}/Reco p_{T}", 200, 0, 50, 500, 0, 2);
  hm->registerHisto(h);

  // reco track w/ truth-track matched vs reco pT histograms with quality cuts
  h = new TH1F(TString(get_histo_prefix()) + "nGen_pTReco_cuts",
               "Gen tracks at reco p_{T} (#geq 2 MVTX, #geq 1 INTT, #geq 20 TPC); Reco p_{T} [GeV/c]", 200, 0.1, 50.5);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "nReco_pTReco_cuts",
               ";Gen p_{T} [GeV/c];Track count / bin", 200, 0.1, 50.5);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "pTRecoTruthMatchedRatio_pTReco_cuts",
               ";Reco p_{T} [GeV/c];Matched p_{T}/Reco p_{T}", 200, 0, 50, 500, 0, 2);
  hm->registerHisto(h);

  // reco pT histogram
  h = new TH1F(TString(get_histo_prefix()) + "nReco_pTGen",
               "Reco tracks at truth p_{T};Truth p_{T} [GeV/c]", 200, 0.1, 50.5);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  // reco pT histogram vs tracker clusters
  h = new TH2F(TString(get_histo_prefix()) + "nMVTX_nReco_pTGen",
               "Reco tracks at truth p_{T};Truth p_{T} [GeV/c];nCluster_{MVTX}", 200, 0.1, 50.5, 6, -.5, 5.5);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "nINTT_nReco_pTGen",
               "Reco tracks at truth p_{T};Truth p_{T} [GeV/c];nCluster_{INTT}", 200, 0.1, 50.5, 6, -.5, 5.5);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "nTPC_nReco_pTGen",
               "Reco tracks at truth p_{T};Truth p_{T} [GeV/c];nCluster_{TPC}", 200, 0.1, 50.5, 60, -.5, 59.5);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  // DCA histograms
  h = new TH2F(TString(get_histo_prefix()) + "DCArPhi_pT",
               "DCA resolution at truth p_{T};Truth p_{T} [GeV/c];DCA(r#phi) resolution [cm]", 200, 0.1, 50.5, 500, -0.05, 0.05);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "DCAZ_pT",
               "DCA resolution at truth p_{T};Truth p_{T} [GeV/c];DCA(Z) resolution [cm]", 200, 0.1, 50.5, 500, -0.05, 0.05);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  // DCA histograms with cuts
  h = new TH2F(TString(get_histo_prefix()) + "DCArPhi_pT_cuts",
               "DCA Resolution (#geq 2 MVTX, #geq 1 INTT, #geq 20 TPC);Truth p_{T} [GeV/c];DCA(r#phi) resolution [cm]", 200, 0.1, 50.5, 500, -0.05, 0.05);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "DCAZ_pT_cuts",
               "DCA Resolution (#geq 2 MVTX, #geq 1 INTT, #geq 20 TPC);Truth p_{T} [GeV/c];DCA(Z) resolution [cm]", 200, 0.1, 50.5, 500, -0.05, 0.05);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "SigmalizedDCArPhi_pT",
               "Sigmalized DCA (#geq 2 MVTX, #geq 1 INTT, #geq 20 TPC);Truth p_{T} [GeV/c];Sigmalized DCA(r#phi) [cm]", 200, 0.1, 50.5, 500, -5., 5.);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH2F(TString(get_histo_prefix()) + "SigmalizedDCAZ_pT",
               "Sigmalized DCA (#geq 2 MVTX, #geq 1 INTT, #geq 20 TPC);Truth p_{T} [GeV/c];Sigmalized DCA(Z) [cm]", 200, 0.1, 50.5, 500, -5., 5.);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  // reco pT histogram
  h = new TH1F(TString(get_histo_prefix()) + "nGen_pTGen",
               ";Truth p_{T} [GeV/c];Track count / bin", 200, 0.1, 50.5);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  // reco pT histogram
  h = new TH1F(TString(get_histo_prefix()) + "nReco_etaGen",
               "Reco tracks at truth  #eta;Truth  #eta", 200, -2, 2);
  //  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  // reco pT histogram
  h = new TH1F(TString(get_histo_prefix()) + "nGen_etaGen",
               ";Truth #eta;Track count / bin", 200, -2, 2);
  //  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  // clusters per layer and per track histogram
  h = new TH1F(TString(get_histo_prefix()) + "nClus_layer", "Reco Clusters per layer per track;Layer;nCluster", 64, 0, 64);
  hm->registerHisto(h);

  // clusters per layer and per generated track histogram
  h = new TH1F(TString(get_histo_prefix()) + "nClus_layerGen", "Reco Clusters per layer per truth track;Layer;nCluster", 64, 0, 64);
  hm->registerHisto(h);

  // n events and n tracks histogram
  h = new TH1F(TString(get_histo_prefix()) + "Normalization",
               TString(get_histo_prefix()) + " Normalization;Items;Count", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event");
  h->GetXaxis()->SetBinLabel(i++, "Truth Track");
  h->GetXaxis()->SetBinLabel(i++, "Truth Track+");
  h->GetXaxis()->SetBinLabel(i++, "Truth Track-");
  h->GetXaxis()->SetBinLabel(i++, "Reco Track");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
}

void QAG4SimulationTracking::addEmbeddingID(int embeddingID)
{
  m_embeddingIDs.insert(embeddingID);
}

int QAG4SimulationTracking::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
  {
    std::cout << "QAG4SimulationTracking::process_event" << std::endl;
  }

  // load relevant nodes from NodeTree
  load_nodes(topNode);

  // histogram manager
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  if (m_svtxEvalStack)
    m_svtxEvalStack->next_event(topNode);

  SvtxTrackEval *trackeval = m_svtxEvalStack->get_track_eval();
  assert(trackeval);
  SvtxTruthEval *trutheval = m_svtxEvalStack->get_truth_eval();
  assert(trutheval);

  // reco pT / gen pT histogram
  TH1 *h_pTRecoGenRatio = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "pTRecoGenRatio"));
  assert(h_pTRecoGenRatio);

  // reco pT / gen pT histogram
  TH2 *h_pTRecoGenRatio_pTGen = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "pTRecoGenRatio_pTGen"));
  assert(h_pTRecoGenRatio);

  // reco track, truth track matched histogram at reco pT
  TH1 *h_nGen_pTReco = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nGen_pTReco"));
  assert(h_nGen_pTReco);
  // Normalization histogram for reco track truth track matched
  TH1 *h_nReco_pTReco = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nReco_pTReco"));
  assert(h_nReco_pTReco);
  // Truth matched ratio histogram
  TH2 *h_pTRecoTruthMatchedRatio_pTReco = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "pTRecoTruthMatchedRatio_pTReco"));
  assert(h_pTRecoTruthMatchedRatio_pTReco);

  // reco track, truth track matched histogram at reco pT with cuts
  TH1 *h_nGen_pTReco_cuts = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nGen_pTReco_cuts"));
  assert(h_nGen_pTReco_cuts);
  // Normalization histogram for reco track truth track matched with cuts
  TH1 *h_nReco_pTReco_cuts = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nReco_pTReco_cuts"));
  assert(h_nReco_pTReco_cuts);
  // Truth matched ratio histogram with cuts
  TH2 *h_pTRecoTruthMatchedRatio_pTReco_cuts = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "pTRecoTruthMatchedRatio_pTReco_cuts"));
  assert(h_pTRecoTruthMatchedRatio_pTReco_cuts);

  // reco histogram plotted at gen pT
  TH1 *h_nReco_pTGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nReco_pTGen"));
  assert(h_nReco_pTGen);

  // reco histogram plotted at gen pT
  TH1 *h_nMVTX_nReco_pTGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nMVTX_nReco_pTGen"));
  assert(h_nMVTX_nReco_pTGen);
  // reco histogram plotted at gen pT
  TH1 *h_nINTT_nReco_pTGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nINTT_nReco_pTGen"));
  assert(h_nINTT_nReco_pTGen);
  // reco histogram plotted at gen pT
  TH1 *h_nTPC_nReco_pTGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nTPC_nReco_pTGen"));
  assert(h_nTPC_nReco_pTGen);

  // DCA resolution histogram
  TH2 *h_DCArPhi_pT = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "DCArPhi_pT"));
  assert(h_DCArPhi_pT);
  // DCA resolution histogram
  TH2 *h_DCAZ_pT = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "DCAZ_pT"));
  assert(h_DCAZ_pT);

  // DCA resolution histogram with cuts
  TH2 *h_DCArPhi_pT_cuts = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "DCArPhi_pT_cuts"));
  assert(h_DCArPhi_pT_cuts);
  // DCA resolution histogram with cuts
  TH2 *h_DCAZ_pT_cuts = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "DCAZ_pT_cuts"));
  assert(h_DCAZ_pT_cuts);
  // Sigmalized DCA histograms with cuts
  TH2 *h_SigmalizedDCArPhi_pT = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "SigmalizedDCArPhi_pT"));
  assert(h_SigmalizedDCArPhi_pT);
  // Sigmalized DCA histograms with cuts
  TH2 *h_SigmalizedDCAZ_pT = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "SigmalizedDCAZ_pT"));
  assert(h_SigmalizedDCAZ_pT);

  // gen pT histogram
  TH1 *h_nGen_pTGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nGen_pTGen"));
  assert(h_nGen_pTGen);

  // reco histogram plotted at gen eta
  TH1 *h_nReco_etaGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nReco_etaGen"));
  assert(h_nReco_etaGen);

  // gen eta histogram
  TH1 *h_nGen_etaGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nGen_etaGen"));
  assert(h_nGen_etaGen);

  // clusters per layer and per track
  auto h_nClus_layer = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nClus_layer"));
  assert(h_nClus_layer);

  // clusters per layer and per generated track
  auto h_nClus_layerGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nClus_layerGen"));
  assert(h_nClus_layer);

  // n events and n tracks histogram
  TH1 *h_norm = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);

  // fill histograms that need truth information
  if (!m_truthContainer)
  {
    std::cout << "QAG4SimulationTracking::process_event - fatal error - missing m_truthContainer! ";
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // build map of clusters associated to G4Particles
  /* inspired from PHTruthTrackSeeding code */
  using KeySet = std::set<TrkrDefs::cluskey>;
  using ParticleMap = std::map<int, KeySet>;
  ParticleMap g4particle_map;

  if (m_cluster_map)
  {
    // loop over clusters
    for (const auto &hitsetkey : m_cluster_map->getHitSetKeys())
    {
      auto range = m_cluster_map->getClusters(hitsetkey);
      for (auto clusterIter = range.first; clusterIter != range.second; ++clusterIter)
      {
        // store cluster key
        const auto &key = clusterIter->first;

        // loop over associated g4hits
        for (const auto &g4hit : find_g4hits(key))
        {
          const int trkid = g4hit->get_trkid();
          auto iter = g4particle_map.lower_bound(trkid);
          if (iter != g4particle_map.end() && iter->first == trkid)
          {
            iter->second.insert(key);
          }
          else
          {
            g4particle_map.insert(iter, std::make_pair(trkid, KeySet({key})));
          }
        }
      }
    }
  }

  // loop over reco tracks to fill norm histogram for track matching
  if (m_trackMap)
  {
    for (SvtxTrackMap::Iter iter = m_trackMap->begin();
         iter != m_trackMap->end();
         ++iter)
    {
      SvtxTrack *track = iter->second;
      if (!track)
      {
        continue;
      }

      const double px = track->get_px();
      const double py = track->get_py();
      const double pz = track->get_pz();
      const TVector3 v(px, py, pz);
      const double pt = v.Pt();
      h_nReco_pTReco->Fill(pt);  // normalization histogram fill

      int MVTX_hits = 0;
      int INTT_hits = 0;
      int TPC_hits = 0;

      TrackSeed *tpcseed = track->get_tpc_seed();
      TrackSeed *silseed = track->get_silicon_seed();
      if (silseed)
      {
        for (auto cluster_iter = silseed->begin_cluster_keys();
             cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
        {
          const auto &cluster_key = *cluster_iter;
          const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

          if (trackerID == TrkrDefs::mvtxId)
            ++MVTX_hits;
          else if (trackerID == TrkrDefs::inttId)
            ++INTT_hits;
        }
      }
      if (tpcseed)
      {
        for (auto cluster_iter = tpcseed->begin_cluster_keys(); cluster_iter != tpcseed->end_cluster_keys(); ++cluster_iter)
        {
          const auto &cluster_key = *cluster_iter;
          const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

          if (trackerID == TrkrDefs::tpcId) ++TPC_hits;
        }
      }

      if (MVTX_hits >= 2 && INTT_hits >= 1 && TPC_hits >= 20)
      {
        h_nReco_pTReco_cuts->Fill(pt);  // normalization histogram fill with cuts
      }

      auto g4particle_match = trackeval->max_truth_particle_by_nclusters(track);
      if (g4particle_match)
      {
        SvtxTrack *matched_track = trackeval->best_track_from(g4particle_match);
        if (matched_track)
        {
          if (matched_track->get_id() == track->get_id())
          {
            h_nGen_pTReco->Fill(pt);  // fill if matching truth track

            const double gpx = g4particle_match->get_px();
            const double gpy = g4particle_match->get_py();
            TVector3 gv(gpx, gpy, 0);
            const double gpt = gv.Pt();

            const double pt_ratio = (pt != 0) ? gpt / pt : 0;
            h_pTRecoTruthMatchedRatio_pTReco->Fill(pt, pt_ratio);

            if (MVTX_hits >= 2 && INTT_hits >= 1 && TPC_hits >= 20)
            {
              h_nGen_pTReco_cuts->Fill(pt);
              h_pTRecoTruthMatchedRatio_pTReco_cuts->Fill(pt, pt_ratio);
            }
          }
        }
      }
    }
  }
  else
  {
    std::cout << __PRETTY_FUNCTION__ << " : Fatal error: missing SvtxTrackMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }  // reco track loop

  PHG4TruthInfoContainer::ConstRange range = m_truthContainer->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    // get the truth particle information
    PHG4Particle *g4particle = iter->second;

    if (Verbosity())
    {
      std::cout << "QAG4SimulationTracking::process_event - processing ";
      g4particle->identify();
    }

    if (!m_embeddingIDs.empty())
    {
      // only analyze subset of particle with proper embedding IDs
      int candidate_embedding_id = trutheval->get_embed(g4particle);
      if (candidate_embedding_id <= m_embed_id_cut) candidate_embedding_id = -1;

      // skip if no match
      if (m_embeddingIDs.find(candidate_embedding_id) == m_embeddingIDs.end()) continue;
    }

    double gpx = g4particle->get_px();
    double gpy = g4particle->get_py();
    double gpz = g4particle->get_pz();
    double gpt = 0;
    double geta = NAN;

    if (gpx != 0 && gpy != 0)
    {
      TVector3 gv(gpx, gpy, gpz);
      gpt = gv.Pt();
      geta = gv.Eta();
      //      gphi = gv.Phi();
    }
    if (m_etaRange.first < geta and geta < m_etaRange.second)
    {
      if (Verbosity())
      {
        std::cout << "QAG4SimulationTracking::process_event - accept particle eta = " << geta << std::endl;
      }
    }
    else
    {
      if (Verbosity())
        std::cout << "QAG4SimulationTracking::process_event - ignore particle eta = " << geta << std::endl;
      continue;
    }

    const int pid = g4particle->get_pid();
    TParticlePDG *pdg_p = TDatabasePDG::Instance()->GetParticle(pid);
    if (!pdg_p)
    {
      std::cout << "QAG4SimulationTracking::process_event - Error - invalid particle ID = " << pid << std::endl;
      continue;
    }

    const double gcharge = pdg_p->Charge() / 3;
    if (gcharge > 0)
    {
      h_norm->Fill("Truth Track+", 1);
    }
    else if (gcharge < 0)
    {
      h_norm->Fill("Truth Track-", 1);
    }
    else
    {
      if (Verbosity())
        std::cout << "QAG4SimulationTracking::process_event - invalid particle ID = " << pid << std::endl;
      continue;
    }
    h_norm->Fill("Truth Track", 1);

    h_nGen_pTGen->Fill(gpt);
    h_nGen_etaGen->Fill(geta);

    // loop over clusters associated to this G4Particle
    {
      const auto mapIter = g4particle_map.find(iter->first);
      if (mapIter != g4particle_map.cend())
      {
        for (const auto &cluster_key : mapIter->second)
        {
          h_nClus_layerGen->Fill(TrkrDefs::getLayer(cluster_key));
        }
      }
      else if (Verbosity())
      {
        std::cout << "QAG4SimulationTracking::process_event - could nof find clusters associated to G4Particle " << iter->first << std::endl;
      }
    }
    // look for best matching track in reco data & get its information
    SvtxTrack *track = trackeval->best_track_from(g4particle);
    if (track)
    {
      bool match_found(false);

      if (m_uniqueTrackingMatch)
      {
        PHG4Particle *g4particle_matched = trackeval->max_truth_particle_by_nclusters(track);
        if (g4particle_matched)
        {
          if (g4particle_matched->get_track_id() == g4particle->get_track_id())
          {
            match_found = true;
            if (Verbosity())
              std::cout << "QAG4SimulationTracking::process_event - found unique match for g4 particle " << g4particle->get_track_id() << std::endl;
          }
          else
          {
            if (Verbosity())
              std::cout << "QAG4SimulationTracking::process_event - none unique match for g4 particle " << g4particle->get_track_id()
                        << ". The track belong to g4 particle " << g4particle_matched->get_track_id() << std::endl;
          }
        }  //        if (g4particle_matched)
        else
        {
          if (Verbosity())
            std::cout << "QAG4SimulationTracking::process_event - none unique match for g4 particle " << g4particle->get_track_id()
                      << ". The track belong to no g4 particle!" << std::endl;
        }
      }

      if (match_found || !m_uniqueTrackingMatch)
      {
        h_nReco_etaGen->Fill(geta);
        h_nReco_pTGen->Fill(gpt);

        float dca3dxy = NAN;
        float dca3dz = NAN;
        float dca3dxysigma = NAN;
        float dca3dzsigma = NAN;
        get_dca(track, dca3dxy, dca3dz, dca3dxysigma, dca3dzsigma);

        double px = track->get_px();
        double py = track->get_py();
        double pz = track->get_pz();
        double pt;
        TVector3 v(px, py, pz);
        pt = v.Pt();
        //        eta = v.Pt();
        //      phi = v.Pt();

        float pt_ratio = (gpt != 0) ? pt / gpt : 0;
        h_pTRecoGenRatio->Fill(pt_ratio);
        h_pTRecoGenRatio_pTGen->Fill(gpt, pt_ratio);
        h_DCArPhi_pT->Fill(pt, dca3dxy);
        h_DCAZ_pT->Fill(pt, dca3dz);
        h_norm->Fill("Reco Track", 1);

        int MVTX_hits = 0;
        int INTT_hits = 0;
        int TPC_hits = 0;

        auto tpcSeed = track->get_tpc_seed();
        auto silSeed = track->get_silicon_seed();

        if (silSeed)
        {
          for (auto cluster_iter = silSeed->begin_cluster_keys();
               cluster_iter != silSeed->end_cluster_keys(); ++cluster_iter)
          {
            const auto &cluster_key = *cluster_iter;
            const auto trackerID = TrkrDefs::getTrkrId(cluster_key);
            const auto layer = TrkrDefs::getLayer(cluster_key);

            h_nClus_layer->Fill(layer);
            if (trackerID == TrkrDefs::mvtxId)
              ++MVTX_hits;
            else if (trackerID == TrkrDefs::inttId)
              ++INTT_hits;
          }
        }

        if (tpcSeed)
        {
          for (auto cluster_iter = tpcSeed->begin_cluster_keys(); cluster_iter != tpcSeed->end_cluster_keys(); ++cluster_iter)
          {
            const auto &cluster_key = *cluster_iter;
            const auto trackerID = TrkrDefs::getTrkrId(cluster_key);
            const auto layer = TrkrDefs::getLayer(cluster_key);

            h_nClus_layer->Fill(layer);
            if (trackerID == TrkrDefs::tpcId) ++TPC_hits;
          }
        }

        if (MVTX_hits >= 2 && INTT_hits >= 1 && TPC_hits >= 20)
        {
          h_DCArPhi_pT_cuts->Fill(pt, dca3dxy);
          h_DCAZ_pT_cuts->Fill(pt, dca3dz);
          h_SigmalizedDCArPhi_pT->Fill(pt, dca3dxy / dca3dxysigma);
          h_SigmalizedDCAZ_pT->Fill(pt, dca3dz / dca3dzsigma);
        }

        h_nMVTX_nReco_pTGen->Fill(gpt, MVTX_hits);
        h_nINTT_nReco_pTGen->Fill(gpt, INTT_hits);
        h_nTPC_nReco_pTGen->Fill(gpt, TPC_hits);
      }  //      if (match_found)

    }  //    if (track)
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void QAG4SimulationTracking::get_dca(SvtxTrack *track, float &dca3dxy,
                                     float &dca3dz, float &dca3dxysigma,
                                     float &dca3dzsigma)
{
  Acts::Vector3 pos(track->get_x(),
                    track->get_y(),
                    track->get_z());
  Acts::Vector3 mom(track->get_px(),
                    track->get_py(),
                    track->get_pz());

  auto vtxid = track->get_vertex_id();
  auto svtxVertex = m_vertexMap->get(vtxid);
  if (!svtxVertex) return;

  Acts::Vector3 vertex(svtxVertex->get_x(),
                       svtxVertex->get_y(),
                       svtxVertex->get_z());

  pos -= vertex;

  Acts::ActsSymMatrix<3> posCov;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      posCov(i, j) = track->get_error(i, j);
    }
  }

  Acts::Vector3 r = mom.cross(Acts::Vector3(0., 0., 1.));
  float phi = atan2(r(1), r(0));

  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;
  rot(0, 0) = cos(phi);
  rot(0, 1) = -sin(phi);
  rot(0, 2) = 0;
  rot(1, 0) = sin(phi);
  rot(1, 1) = cos(phi);
  rot(1, 2) = 0;
  rot(2, 0) = 0;
  rot(2, 1) = 0;
  rot(2, 2) = 1;

  rot_T = rot.transpose();

  Acts::Vector3 pos_R = rot * pos;
  Acts::ActsSymMatrix<3> rotCov = rot * posCov * rot_T;

  dca3dxy = pos_R(0);
  dca3dz = pos_R(2);
  dca3dxysigma = sqrt(rotCov(0, 0));
  dca3dzsigma = sqrt(rotCov(2, 2));

  return;
}

int QAG4SimulationTracking::load_nodes(PHCompositeNode *topNode)
{
  m_truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthContainer)
  {
    std::cout << "QAG4SimulationTracking::load_nodes - Fatal Error - "
              << "unable to find DST node "
              << "G4TruthInfo" << std::endl;
    assert(m_truthContainer);
  }

  // cluster map
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  //  assert(m_cluster_map);

  // cluster hit association map
  m_cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  //  assert(m_cluster_hit_map);

  // cluster hit association map
  m_hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  //  assert(m_hit_truth_map);

  // g4hits
  m_g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  m_g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  m_g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  m_g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string
QAG4SimulationTracking::get_histo_prefix()
{
  return std::string("h_") + Name() + std::string("_");
}

QAG4SimulationTracking::G4HitSet QAG4SimulationTracking::find_g4hits(TrkrDefs::cluskey cluster_key) const
{
  assert(m_cluster_hit_map);

  // find hitset associated to cluster
  G4HitSet out;
  const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);

  // loop over hits associated to clusters
  const auto range = m_cluster_hit_map->getHits(cluster_key);
  for (auto iter = range.first; iter != range.second; ++iter)
  {
    // hit key
    const auto &hit_key = iter->second;

    // store hits to g4hit associations
    TrkrHitTruthAssoc::MMap g4hit_map;
    m_hit_truth_map->getG4Hits(hitset_key, hit_key, g4hit_map);

    // find corresponding g4 hist
    for (auto truth_iter = g4hit_map.begin(); truth_iter != g4hit_map.end(); ++truth_iter)
    {
      const auto g4hit_key = truth_iter->second.second;
      PHG4Hit *g4hit = nullptr;

      switch (TrkrDefs::getTrkrId(hitset_key))
      {
      case TrkrDefs::mvtxId:
        assert(m_g4hits_mvtx);
        if (m_g4hits_mvtx) g4hit = m_g4hits_mvtx->findHit(g4hit_key);
        break;

      case TrkrDefs::inttId:
        assert(m_g4hits_intt);
        if (m_g4hits_intt) g4hit = m_g4hits_intt->findHit(g4hit_key);
        break;

      case TrkrDefs::tpcId:
        assert(m_g4hits_tpc);
        if (m_g4hits_tpc) g4hit = m_g4hits_tpc->findHit(g4hit_key);
        break;

      case TrkrDefs::micromegasId:
        assert(m_g4hits_micromegas);
        if (m_g4hits_micromegas) g4hit = m_g4hits_micromegas->findHit(g4hit_key);
        break;

      default:
        break;
      }

      if (g4hit) out.insert(g4hit);
    }
  }

  return out;
}
