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

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <TAxis.h>
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
  hm->registerHisto(new TH1F(TString(get_histo_prefix()) + "pTRecoGenRatio",
                             "Reco p_{T}/Truth p_{T}",500, 0, 2));
  // reco pT histogram
  hm->registerHisto(new TH1F(TString(get_histo_prefix()) + "nReco_pTGen",
                             "Reco tracks at truth p_{T}",200, -0.5, 50.5));
  // reco pT histogram
  hm->registerHisto(new TH1F(TString(get_histo_prefix()) + "nGen_pTGen",
                             "Truth p_{T}",200, -0.5, 50.5));
  
  // n events and n tracks histogram
  TH1F *h = new TH1F(TString(get_histo_prefix()) + "Normalization",
                     TString(get_histo_prefix()) + " Normalization;Items;Count", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event");
  h->GetXaxis()->SetBinLabel(i++, "Track");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
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

  // reco pT / gen pT histogram
  TH1F *h_pTRecoGenRatio = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "pTRecoGenRatio"));
  assert(h_pTRecoGenRatio);
  
  // reco histogram plotted at gen pT
  TH1F *h_nReco_pTGen = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "nReco_pTGen"));
  assert(h_nReco_pTGen);
  
  // gen pT histogram
  TH1F *h_nGen_pTGen = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "nGen_pTGen"));
  assert(h_nGen_pTGen);
  
  // n events and n tracks histogram
  TH1F *h_norm = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);
  
  // fill histograms that need truth information
  if(_truthContainer)
  {
    PHG4TruthInfoContainer::ConstRange range = _truthContainer->GetPrimaryParticleRange();
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
    {
      // get the truth particle information
      PHG4Particle* g4particle = iter->second;
      
      // is this needed? what does it mean?
//      if (_scan_for_embedded)
//      {
//        if (trutheval->get_embed(g4particle) <= 0) continue;
//      }
      
      float gpx = g4particle->get_px();
      float gpy = g4particle->get_py();
      float gpz = g4particle->get_px();
      float gpt = 0;

      if (gpx != 0 && gpy != 0) {
        TVector3 gv(gpx,gpy,gpz);
        gpt = gv.Pt();
  //      geta = gv.Pt();
  //      gphi = gv.Pt();
      }
      h_nGen_pTGen->Fill(gpt);

      // look for best matching track in reco data & get its information
      SvtxTrack* track = trackeval->best_track_from(g4particle);
      if(track) {
        h_nReco_pTGen->Fill(gpt);
        
        float px = track->get_px();
        float py = track->get_py();
        float pz = track->get_pz();
        float pt;
        TVector3 v(px,py,pz);
        pt = v.Pt();
  //      eta = v.Pt();
  //      phi = v.Pt();
        
        float pt_ratio = (gpt != 0) ? pt/gpt : 0;
        h_pTRecoGenRatio->Fill(pt_ratio);
        h_norm->Fill("Track", 1);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

string
QAG4SimulationTracking::get_histo_prefix()
{
  return "h_QAG4Sim_Tracking_";
}
