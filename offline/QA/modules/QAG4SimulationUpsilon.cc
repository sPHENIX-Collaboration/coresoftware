#include "QAG4SimulationUpsilon.h"

#include "QAHistManagerDef.h"

#include <g4eval/SvtxEvalStack.h>

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

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
#include <TParticlePDG.h>  // for TParticlePDG
#include <TString.h>
#include <TVector3.h>

#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/ThreeVector.h>  // for Hep3Vector

#include <cassert>
#include <cmath>
#include <cstdlib>  // for abs
#include <iostream>
#include <utility>  // for pair

QAG4SimulationUpsilon::QAG4SimulationUpsilon(const std::string &name)
  : SubsysReco(name)
  , _svtxEvalStack(nullptr)
  , m_etaRange(-1, 1)
  , _truthContainer(nullptr)
{
}

int QAG4SimulationUpsilon::InitRun(PHCompositeNode *topNode)
{
  _truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                               "G4TruthInfo");
  if (!_truthContainer)
  {
    std::cout << "QAG4SimulationUpsilon::InitRun - Fatal Error - "
              << "unable to find DST node "
              << "G4TruthInfo" << std::endl;
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

int QAG4SimulationUpsilon::Init(PHCompositeNode * /*topNode*/)
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

  h = new TH1F(TString(get_histo_prefix()) + "nGen_Pair_InvMassGen",
               ";Truth Invariant Mass [GeV/c^2];Pair count / bin", 450, 0, 15);
  //  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "nReco_Pair_InvMassReco",
               ";Reco Invariant Mass [GeV/c^2];Pair count / bin", 450, 0, 15);
  //  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);

  // n events and n tracks histogram
  h = new TH1F(TString(get_histo_prefix()) + "Normalization",
               TString(get_histo_prefix()) + " Normalization;Items;Count", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event");
  h->GetXaxis()->SetBinLabel(i++, "Truth Track+");
  h->GetXaxis()->SetBinLabel(i++, "Truth Track-");
  h->GetXaxis()->SetBinLabel(i++, "Reco Track");
  h->GetXaxis()->SetBinLabel(i++, "Truth Upsilon");
  h->GetXaxis()->SetBinLabel(i++, "Truth Upsilon in Acc.");
  h->GetXaxis()->SetBinLabel(i++, "Reco Upsilon");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
}

void QAG4SimulationUpsilon::addEmbeddingID(int embeddingID)
{
  m_embeddingIDs.insert(embeddingID);
}

int QAG4SimulationUpsilon::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
    std::cout << "QAG4SimulationUpsilon::process_event() entered" << std::endl;

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

  // inv mass
  TH1 *h_nGen_Pair_InvMassGen = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nGen_Pair_InvMassGen"));
  assert(h_nGen_etaGen);
  // inv mass
  TH1 *h_nReco_Pair_InvMassReco = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "nReco_Pair_InvMassReco"));
  assert(h_nGen_etaGen);

  // n events and n tracks histogram
  TH1 *h_norm = dynamic_cast<TH1 *>(hm->getHisto(get_histo_prefix() + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);

  // buffer for daugther particles
  typedef std::set<std::pair<PHG4Particle *, SvtxTrack *>> truth_reco_set_t;
  truth_reco_set_t truth_reco_set_pos;
  truth_reco_set_t truth_reco_set_neg;

  // fill histograms that need truth information
  if (!_truthContainer)
  {
    std::cout << "QAG4SimulationUpsilon::process_event - fatal error - missing _truthContainer! ";
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHG4TruthInfoContainer::ConstRange range = _truthContainer->GetPrimaryParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    // get the truth particle information
    PHG4Particle *g4particle = iter->second;

    if (Verbosity())
    {
      std::cout << "QAG4SimulationUpsilon::process_event - processing ";
      g4particle->identify();
    }

    if (m_embeddingIDs.size() > 0)
    {
      // only analyze subset of particle with proper embedding IDs
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
        std::cout << "QAG4SimulationUpsilon::process_event - accept particle eta = " << geta << std::endl;
      }
    }
    else
    {
      if (Verbosity())
        std::cout << "QAG4SimulationUpsilon::process_event - ignore particle eta = " << geta << std::endl;
      continue;
    }

    const int pid = g4particle->get_pid();

    if (abs(pid) != abs(m_daughterAbsPID))
    {
      if (Verbosity())
        std::cout << "QAG4SimulationUpsilon::process_event - ignore particle PID = " << pid << " as m_daughterAbsPID = " << m_daughterAbsPID << std::endl;
      continue;
    }

    TParticlePDG *pdg_p = TDatabasePDG::Instance()->GetParticle(pid);
    if (!pdg_p)
    {
      std::cout << "QAG4SimulationUpsilon::process_event - Error - invalid particle ID = " << pid << std::endl;
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
        std::cout << "QAG4SimulationUpsilon::process_event - invalid neutral decay particle ID = " << pid << std::endl;
      continue;
    }
    //        h_norm->Fill("Truth Track", 1);

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

    if (gcharge > 0)
    {
      truth_reco_set_pos.insert(std::make_pair(g4particle, track));
    }
    else if (gcharge < 0)
    {
      truth_reco_set_neg.insert(std::make_pair(g4particle, track));
    }
  }  //  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)

  // building pairs with buffer for daugter particles
  TParticlePDG *pdg_p = TDatabasePDG::Instance()->GetParticle(m_quarkoniaPID);
  if (!pdg_p)
  {
    std::cout << "QAG4SimulationUpsilon::process_event - Fatal Error - invalid particle ID m_quarkoniaPID = " << m_quarkoniaPID << std::endl;

    return Fun4AllReturnCodes::ABORTRUN;
  }
  const double quarkonium_mass = pdg_p->Mass();
  TParticlePDG *pdg_d = TDatabasePDG::Instance()->GetParticle(m_daughterAbsPID);
  if (!pdg_d)
  {
    std::cout << "QAG4SimulationUpsilon::process_event - Fatal Error - invalid particle ID m_daughterAbsPID = " << m_daughterAbsPID << std::endl;

    return Fun4AllReturnCodes::ABORTRUN;
  }
  const double daughter_mass = pdg_d->Mass();

  for (const auto &pair_pos : truth_reco_set_pos)
    for (const auto &pair_neg : truth_reco_set_neg)
    {
      assert(pair_pos.first);
      assert(pair_neg.first);

      const CLHEP::HepLorentzVector gv_pos(
          pair_pos.first->get_px(),
          pair_pos.first->get_py(),
          pair_pos.first->get_pz(),
          pair_pos.first->get_e());

      const CLHEP::HepLorentzVector gv_neg(
          pair_neg.first->get_px(),
          pair_neg.first->get_py(),
          pair_neg.first->get_pz(),
          pair_neg.first->get_e());

      const CLHEP::HepLorentzVector gv_quakonium = gv_pos + gv_neg;

      if (fabs(quarkonium_mass - gv_quakonium.m()) > 1e-3)
      {
        if (Verbosity())
        {
          std::cout << "QAG4SimulationUpsilon::process_event - invalid pair with in compativle mass with " << quarkonium_mass << "GeV for PID = " << m_quarkoniaPID << ": " << std::endl;
          pair_pos.first->identify();
          pair_neg.first->identify();
        }
        continue;
      }

      h_nGen_Pair_InvMassGen->Fill(gv_quakonium.m());
      h_norm->Fill("Truth Upsilon in Acc.", 1);

      if (pair_pos.second and pair_neg.second)
      {
        CLHEP::HepLorentzVector v_pos;
        CLHEP::HepLorentzVector v_neg;

        v_pos.setVectM(
            CLHEP::Hep3Vector(
                pair_pos.second->get_px(),
                pair_pos.second->get_py(),
                pair_pos.second->get_pz()),
            daughter_mass);
        v_neg.setVectM(
            CLHEP::Hep3Vector(
                pair_neg.second->get_px(),
                pair_neg.second->get_py(),
                pair_neg.second->get_pz()),
            daughter_mass);

        const CLHEP::HepLorentzVector v_quakonium = v_pos + v_neg;

        h_nReco_Pair_InvMassReco->Fill(v_quakonium.m());
        h_norm->Fill("Reco Upsilon", 1);

      }  //      if (pair_pos.second and pair_neg.second)

    }  //    for (const auto &pair_neg : truth_reco_set_neg)

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string
QAG4SimulationUpsilon::get_histo_prefix()
{
  return std::string("h_") + Name() + std::string("_");
}
