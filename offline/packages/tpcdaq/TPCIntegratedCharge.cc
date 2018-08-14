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

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <CLHEP/Units/SystemOfUnits.h>

#include <boost/format.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;
using namespace CLHEP;

TPCIntegratedCharge::TPCIntegratedCharge(
    unsigned int minLayer,
    unsigned int m_maxLayer,
    const std::string& outputfilename)
  : SubsysReco("TPCIntegratedCharge")
  , m_outputFileName(outputfilename)
  , m_minLayer(minLayer)
  , m_maxLayer(m_maxLayer)
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
  //  PHG4CylinderCellGeomContainer* seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  //  if (!seggeo)
  //  {
  //    cout << "could not locate geo node " << "CYLINDERCELLGEOM_SVTX" << endl;
  //    exit(1);
  //  }

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
  h->GetXaxis()->SetBinLabel(i++, "TPC Charge e");
  h->GetXaxis()->SetBinLabel(i++, "TPC Charge fC");
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  //  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  //  {
  //    const PHG4CylinderCellGeom* layer_geom = seggeo->GetLayerCellGeom(layer);

  //    const string histNameCellHit(boost::str(boost::format{"hCellHit_Layer%1%"} % layer));
  //    const string histNameCellCharge(boost::str(boost::format{"hCellCharge_Layer%1%"} % layer));

  //  }

  hm->registerHisto(new TH2D("hLayerCellHit",  //
                             "Number of ADC time-bin hit per channel;Layer ID;Hit number",
                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                             300, -.5, 299.5));
  hm->registerHisto(new TH2D("hLayerCellCharge",  //
                             "Charge integrated over drift window per channel;Layer ID;Charge [fC]",
                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                             1000, 0, 1e7 * eplus / (1e-15 * coulomb)));

  hm->registerHisto(new TH2D("hLayerSumCellHit",  //
                             "Number of ADC time-bin hit integrated over channels per layer;Layer ID;Hit number",
                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                             10000, -.5, 99999.5));
  hm->registerHisto(new TH2D("hLayerSumCellCharge",  //
                             "Charge integrated over drift window and channel per layer;Layer ID;Charge [fC]",
                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                             10000, 0, 1000 * 4e6 * eplus / (1e-15 * coulomb)));

  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCIntegratedCharge::process_event(PHCompositeNode* topNode)
{
  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);
  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto("hNormalization"));
  assert(h_norm);

  PHG4HitContainer* g4hit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_SVTX");
  if (!g4hit)
  {
    cout << "TPCIntegratedCharge::process_event - Could not locate g4 hit node G4HIT_SVTX" << endl;
    exit(1);
  }
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_SVTX");
  if (!cells)
  {
    cout << "TPCIntegratedCharge::process_event - could not locate cell node "
         << "G4CELL_SVTX" << endl;
    exit(1);
  }

  PHG4CylinderCellGeomContainer* seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!seggeo)
  {
    cout << "TPCIntegratedCharge::process_event - could not locate geo node "
         << "CYLINDERCELLGEOM_SVTX" << endl;
    exit(1);
  }
  PHHepMCGenEventMap* geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    static bool once = true;
    if (once)
    {
      once = false;

      std::cout << "TPCIntegratedCharge::process_event - missing node PHHepMCGenEventMap. Skipping HepMC stat." << std::endl;
    }
  }
  else
  {
    h_norm->Fill("Collision count", geneventmap->size());
  }

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

  // prepreare charge stat.
  unsigned int nZBins = 0;
  vector<array<vector<int>, 2> > layerChanCellHit(m_maxLayer + 1);
  vector<array<vector<double>, 2> > layerChanCellCharge(m_maxLayer + 1);
  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    PHG4CylinderCellGeom* layerGeom =
        seggeo->GetLayerCellGeom(layer);
    assert(layerGeom);

    // start with an empty vector of cells for each phibin
    const int nphibins = layerGeom->get_phibins();
    assert(nphibins > 0);

    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "TPCIntegratedCharge::process_event - init layer " << layer << " with "
           << "nphibins = " << nphibins
           << ", layerGeom->get_zbins() = " << layerGeom->get_zbins() << endl;
    }

    if (nZBins <= 0)
    {
      nZBins = layerGeom->get_zbins();
      assert(nZBins > 0);
    }
    else
    {
      if ((int) nZBins != layerGeom->get_zbins())
      {
        cout << "TPCIntegratedCharge::process_event - Fatal Error - nZBin at layer " << layer << " is " << layerGeom->get_zbins()
             << ", which is different from previous layers of nZBin = " << nZBins << endl;
        exit(1);
      }
    }

    layerChanCellHit[layer][0].resize(nphibins, 0);
    layerChanCellHit[layer][1].resize(nphibins, 0);

    layerChanCellCharge[layer][0].resize(nphibins, 0);
    layerChanCellCharge[layer][1].resize(nphibins, 0);
  }
  assert(nZBins > 0);

  // count all cells
  PHG4CellContainer::ConstRange cellrange = cells->getCells();
  for (PHG4CellContainer::ConstIterator celliter = cellrange.first;
       celliter != cellrange.second;
       ++celliter)
  {
    PHG4Cell* cell = celliter->second;
    const unsigned int layer = cell->get_layer();
    const unsigned int phibin = PHG4CellDefs::SizeBinning::get_phibin(cell->get_cellid());
    const unsigned int zbin = PHG4CellDefs::SizeBinning::get_zbin(cell->get_cellid());

    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "TPCIntegratedCharge::process_event - process cell: "
           << "layer = " << layer
           << ", get_phibin = " << phibin
           << ", get_zbin = " << zbin
           << ", get_edep = " << cell->get_edep()
           << ", get_eion = " << cell->get_eion() << ": "
           << "Class = " << cell->ClassName();
      cell->identify();
    }

    if (layer < m_minLayer or layer > m_maxLayer) continue;

    if (not(zbin >= 0))
    {
      cout << "TPCIntegratedCharge::process_event - invalid cell: "
           << "layer = " << layer
           << ", get_phibin = " << phibin
           << ", get_zbin = " << zbin
           << ", get_edep = " << cell->get_edep()
           << ", get_eion = " << cell->get_eion() << ": "
           << "Class = " << cell->ClassName();
      cell->identify();
    }

    assert(zbin >= 0);
    assert(zbin < nZBins);
    const int side = (zbin < nZBins / 2) ? 0 : 1;

    vector<int>& ChanCellHit = layerChanCellHit[layer][side];
    vector<double>& ChanCellCharge = layerChanCellCharge[layer][side];

    assert(phibin >= 0);
    assert(phibin < ChanCellHit.size());
    assert(phibin < ChanCellHit.size());

    ChanCellHit[phibin] += 1;
    ChanCellCharge[phibin] += cell->get_edep();
  }

  // fill histograms
  TH2* hLayerCellHit = dynamic_cast<TH2*>(hm->getHisto("hLayerCellHit"));
  assert(hLayerCellHit);
  TH2* hLayerCellCharge = dynamic_cast<TH2*>(hm->getHisto("hLayerCellCharge"));
  assert(hLayerCellCharge);
  TH2* hLayerSumCellHit = dynamic_cast<TH2*>(hm->getHisto("hLayerSumCellHit"));
  assert(hLayerSumCellHit);
  TH2* hLayerSumCellCharge = dynamic_cast<TH2*>(hm->getHisto("hLayerSumCellCharge"));
  assert(hLayerSumCellCharge);

  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    for (unsigned int side = 0; side < 2; ++side)
    {
      int sumHit = 0;
      for (unsigned int phibin = 0; phibin < layerChanCellHit[layer][side].size(); ++phibin)
      {
        const int& hit = layerChanCellHit[layer][side][phibin];
        sumHit += hit;

        hLayerCellHit->Fill(layer, hit);
        h_norm->Fill("TPC Pad Hit", hit);

        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << "TPCIntegratedCharge::process_event - hit: "
               << "hit = " << hit
               << "at layer = " << layer
               << ", side = " << side
               << ", phibin = " << phibin
               << endl;
        }
      }

      if (Verbosity() >= VERBOSITY_MORE)
      {
        cout << "TPCIntegratedCharge::process_event - hLayerSumCellHit->Fill(" << layer << ", " << sumHit << ")" << endl;
      }

      hLayerSumCellHit->Fill(layer, sumHit);

      double sum_charge_fC = 0;
      for (unsigned int phibin = 0; phibin < layerChanCellCharge[layer][side].size(); ++phibin)
      {
        const double& charge_e = layerChanCellCharge[layer][side][phibin];
        const double charge_fC = charge_e * eplus / (1e-15 * coulomb);
        sum_charge_fC += charge_fC;

        hLayerCellCharge->Fill(layer, charge_fC);

        h_norm->Fill("TPC Charge e", charge_e);
        h_norm->Fill("TPC Charge fC", charge_fC);
      }
      hLayerSumCellCharge->Fill(layer, sum_charge_fC);
    }  //    for (unsigned int side = 0; side < 2; ++side)

  }  //  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  h_norm->Fill("Event count", 1);

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
