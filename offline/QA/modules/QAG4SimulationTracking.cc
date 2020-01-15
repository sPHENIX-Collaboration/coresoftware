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

#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>  // for reverse_iterator
#include <utility>   // for pair
#include <vector>

using namespace std;

QAG4SimulationTracking::QAG4SimulationTracking()
  : SubsysReco("QAG4SimulationTracking")
  , m_svtxEvalStack(nullptr)
  , m_truthContainer(nullptr)
{
}

int QAG4SimulationTracking::InitRun(PHCompositeNode *topNode)
{
  m_truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                "G4TruthInfo");
  if (!m_truthContainer)
  {
    cout << "QAG4SimulationTracking::InitRun - Fatal Error - "
         << "unable to find DST node "
         << "G4TruthInfo" << endl;
    assert(m_truthContainer);
  }

  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack.reset(new SvtxEvalStack(topNode));
    m_svtxEvalStack->set_strict(true);
    m_svtxEvalStack->set_verbosity(Verbosity() + 1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationTracking::Init(PHCompositeNode *topNode)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D *h = new TH1D(TString(get_histo_prefix()) + "Normalization",  //
                     TString(get_histo_prefix()) + " Normalization;Items;Count", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event");
  h->GetXaxis()->SetBinLabel(i++, "Track");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  hm->registerHisto(
      new TH1F(
          TString(get_histo_prefix()) + "_pT_gpTReco",  //
          "Reco p_{T}/Truth p_{T}",
          1000, 0, 2));

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationTracking::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
    cout << "QAG4SimulationTracking::process_event() entered" << endl;

  // histogram manager
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  if (m_svtxEvalStack)
    m_svtxEvalStack->next_event(topNode);

  // fill histogram
  //  PHG4Particle *primary = get_truth_particle();
  //  if (!primary)
  //    return Fun4AllReturnCodes::DISCARDEVENT;

  SvtxTrackEval *trackeval = m_svtxEvalStack->get_track_eval();
  assert(trackeval);

  TH1F *h_pT_gpTReco = dynamic_cast<TH1F *>(hm->getHisto(
      get_histo_prefix() + "_pT_gpTReco"  //
      ));
  assert(h_pT_gpTReco);

  h_pT_gpTReco->Fill(1);

  TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto(
      get_histo_prefix() + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);

  return Fun4AllReturnCodes::EVENT_OK;
}

string
QAG4SimulationTracking::get_histo_prefix()
{
  return "h_QAG4Sim_Tracking_";
}
