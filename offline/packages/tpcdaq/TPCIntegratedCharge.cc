#include "TPCIntegratedCharge.h"
#include "TPCDaqDefs.h"

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;

TPCIntegratedCharge::TPCIntegratedCharge(const std::string& outputfilename)
  : SubsysReco("TPCIntegratedCharge")
  , m_outputFileName(outputfilename)
  , m_minLayer(7)
  , m_maxLayer(7 + 48)
{
}

TPCIntegratedCharge::~TPCIntegratedCharge()
{
}

int TPCIntegratedCharge::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCIntegratedCharge::End(PHCompositeNode* topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
    cout << "TPCIntegratedCharge::End - write to " << m_outputFileName << endl;
  PHTFileServer::get().cd(m_outputFileName);

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);
  for (unsigned int i = 0; i < hm->nHistos(); i++)
    hm->getHisto(i)->Write();

  // help index files with TChain
  TTree* T_Index = new TTree("T_Index", "T_Index");
  assert(T_Index);
  T_Index->Write();

  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCIntegratedCharge::InitRun(PHCompositeNode* topNode)
{
  PHG4CylinderCellGeomContainer* seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!seggeo)
  {
    cout << "could not locate geo node " << "CYLINDERCELLGEOM_SVTX" << endl;
    exit(1);
  }

  if (Verbosity() >= VERBOSITY_SOME)
    cout << "TPCIntegratedCharge::get_HistoManager - Making PHTFileServer " << m_outputFileName
         << endl;
  PHTFileServer::get().open(m_outputFileName, "RECREATE");

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);

  TH1D* h = new TH1D("hNormalization",  //
                     "Normalization;Items;Summed quantity", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event count");
  h->GetXaxis()->SetBinLabel(i++, "Collision count");
  h->GetXaxis()->SetBinLabel(i++, "TPC G4Hit");
  h->GetXaxis()->SetBinLabel(i++, "TPC G4Hit Edep");
  h->GetXaxis()->SetBinLabel(i++, "TPC Pad Hit");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCIntegratedCharge::process_event(PHCompositeNode* topNode)
{
  PHG4HitContainer* g4hit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_SVTX");
  if (!g4hit)
  {
    cout << "Could not locate g4 hit node G4HIT_SVTX" << endl;
    exit(1);
  }
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_SVTX");
  if (!cells)
  {
    cout << "could not locate cell node " << "G4CELL_SVTX" << endl;
    exit(1);
  }

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);
  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto("hNormalization"));
  assert(h_norm);

  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(layer);

    for (PHG4HitContainer::ConstIterator hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      const double edep = hiter->second->get_edep();
      h_norm->Fill("TPC G4Hit Edep", edep);
      h_norm->Fill("TPC G4Hit", 1);
    }  //     for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)

  }  //   for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  // at the end, count success events
  h_norm->Fill("Event", 1);

  return Fun4AllReturnCodes::EVENT_OK;
}

Fun4AllHistoManager*
TPCIntegratedCharge::getHistoManager()
{
  static string histname("TPCIntegratedCharge_HISTOS");

  Fun4AllServer* se = Fun4AllServer::instance();
  Fun4AllHistoManager* hm = se->getHistoManager(histname);

  if (not hm)
  {
    cout
        << "TPCIntegratedCharge::get_HistoManager - Making Fun4AllHistoManager " << histname
        << endl;
    hm = new Fun4AllHistoManager(histname);
    se->registerHistoManager(hm);
  }

  assert(hm);

  return hm;
}
