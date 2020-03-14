#include "QAG4SimulationTracking.h"
#include "QAHistManagerDef.h"

#include <g4eval/CaloEvalStack.h>
#include <g4eval/CaloRawClusterEval.h>
#include <g4eval/SvtxEvalStack.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <calobase/RawCluster.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <trackbase_historic/SvtxTrack.h>

#include <g4eval/SvtxTrackEval.h>  // for SvtxTrackEval
#include <g4eval/SvtxTruthEval.h>  // for SvtxTruthEval

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <TAxis.h>
#include <TDatabasePDG.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>
#include <TVector3.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>  // for reverse_iterator
#include <utility>   // for pair
#include <vector>

using namespace std;

QAG4SimulationTracking::QAG4SimulationTracking()
  : SubsysReco("QAG4SimulationTracking")
  , _svtxEvalStack(nullptr)
  , m_etaRange(-1, 1)
  , _truthContainer(nullptr)
{
}

int QAG4SimulationTracking::InitRun(PHCompositeNode *topNode)
{
  _truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                               "G4TruthInfo");
  if (!_truthContainer)
  {
    cout << "QAG4SimulationTracking::InitRun - Fatal Error - "
         << "unable to find DST node "
         << "G4TruthInfo" << endl;
    assert(_truthContainer);
  }

  if (!_svtxEvalStack)
  {
    _svtxEvalStack.reset(new SvtxEvalStack(topNode));
    _svtxEvalStack->set_strict(true);
    _svtxEvalStack->set_verbosity(Verbosity() + 1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationTracking::Init(PHCompositeNode *topNode)
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

  // reco pT histogram
  h = new TH1F(TString(get_histo_prefix()) + "nReco_pTGen",
               "Reco tracks at truth p_{T};Truth p_{T} [GeV/c]", 200, 0.1, 50.5);
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
    cout << "QAG4SimulationTracking::process_event() entered" << endl;

  // histogram manager
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  if (_svtxEvalStack)
    _svtxEvalStack->next_event(topNode);

  SvtxTrackEval *trackeval = _svtxEvalStack->get_track_eval();
  assert(trackeval);
  SvtxTruthEval *trutheval = _svtxEvalStack->get_truth_eval();
  assert(trutheval);

  // reco pT / gen pT histogram
  TH1 *h_pTRecoGenRatio = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "pTRecoGenRatio"));
  assert(h_pTRecoGenRatio);

  // reco pT / gen pT histogram
  TH2 *h_pTRecoGenRatio_pTGen = dynamic_cast<TH2 *>(hm->getHisto(get_histo_prefix() + "pTRecoGenRatio_pTGen"));
  assert(h_pTRecoGenRatio);

  // reco histogram plotted at gen pT
  TH1 *h_nReco_pTGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nReco_pTGen"));
  assert(h_nReco_pTGen);

  // gen pT histogram
  TH1 *h_nGen_pTGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nGen_pTGen"));
  assert(h_nGen_pTGen);

  // reco histogram plotted at gen eta
  TH1 *h_nReco_etaGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nReco_etaGen"));
  assert(h_nReco_etaGen);

  // gen eta histogram
  TH1 *h_nGen_etaGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nGen_etaGen"));
  assert(h_nGen_etaGen);

  // n events and n tracks histogram
  TH1 *h_norm = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);

  // fill histograms that need truth information
  if (_truthContainer)
  {
    PHG4TruthInfoContainer::ConstRange range = _truthContainer->GetPrimaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
    {
      // get the truth particle information
      PHG4Particle *g4particle = iter->second;

      if (Verbosity())
      {
        cout << "QAG4SimulationTracking::process_event - processing ";
        g4particle->identify();
      }

      if (m_embeddingIDs.size() > 0)
      {
        //only analyze subset of particle with proper embedding IDs
        int candidate_embedding_id = trutheval->get_embed(g4particle);
        if (candidate_embedding_id < 0) candidate_embedding_id = -1;

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
          cout << "QAG4SimulationTracking::process_event - accept particle eta = " << geta << endl;
        }
      }
      else
      {
        if (Verbosity())
          cout << "QAG4SimulationTracking::process_event - ignore particle eta = " << geta << endl;
        continue;
      }

      const int pid = g4particle->get_pid();
      TParticlePDG *pdg_p = TDatabasePDG::Instance()->GetParticle(pid);
      if (!pdg_p)
      {
        cout << "QAG4SimulationTracking::process_event - Error - invalid particle ID = " << pid << endl;
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
          cout << "QAG4SimulationTracking::process_event - invalid particle ID = " << pid << endl;
        continue;
      }
      h_norm->Fill("Truth Track", 1);

      h_nGen_pTGen->Fill(gpt);
      h_nGen_etaGen->Fill(geta);

      // look for best matching track in reco data & get its information
      SvtxTrack *track = trackeval->best_track_from(g4particle);
      if (track)
      {
        h_nReco_etaGen->Fill(geta);
        h_nReco_pTGen->Fill(gpt);

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
        h_norm->Fill("Reco Track", 1);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

string
QAG4SimulationTracking::get_histo_prefix()
{
  return string("h_") + Name() + string("_");
}
