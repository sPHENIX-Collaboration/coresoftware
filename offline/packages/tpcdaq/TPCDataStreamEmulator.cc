// $Id: $

/*!
 * \file TPCDataStreamEmulator.cpp
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "TPCDataStreamEmulator.h"

#include "TPCDaqDefs.h"

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <trackbase_historic/SvtxHit.h>
#include <trackbase_historic/SvtxHitMap.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <HepMC/GenEvent.h>

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

TPCDataStreamEmulator::TPCDataStreamEmulator(
    unsigned int minLayer,
    unsigned int m_maxLayer,
    const std::string& outputfilename)
  : SubsysReco("TPCDataStreamEmulator")
  , m_saveDataStreamFile(true)
  , m_outputFileNameBase(outputfilename)
  , m_minLayer(minLayer)
  , m_maxLayer(m_maxLayer)
  , m_evtCounter(-1)
  , m_vertexZAcceptanceCut(10)
  , m_etaAcceptanceCut(1.1)
  , m_hDataSize(nullptr)
  , m_hWavelet(nullptr)
  , m_hNChEta(nullptr)
  , m_hLayerWaveletSize(nullptr)
  , m_hLayerHit(nullptr)
  , m_hLayerZBinHit(nullptr)
  , m_hLayerZBinADC(nullptr)
  , m_hLayerDataSize(nullptr)
  , m_hLayerSumHit(nullptr)
  , m_hLayerSumDataSize(nullptr)
{
}

TPCDataStreamEmulator::~TPCDataStreamEmulator()
{
}

int TPCDataStreamEmulator::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCDataStreamEmulator::End(PHCompositeNode* topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
    cout << "TPCDataStreamEmulator::End - write to " << m_outputFileNameBase + ".root" << endl;
  PHTFileServer::get().cd(m_outputFileNameBase + ".root");

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

int TPCDataStreamEmulator::InitRun(PHCompositeNode* topNode)
{
  PHG4CylinderCellGeomContainer* seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!seggeo)
  {
    cout << "could not locate geo node "
         << "CYLINDERCELLGEOM_SVTX" << endl;
    exit(1);
  }

  int nZBins = 0;
  for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    PHG4CylinderCellGeom* layerGeom =
        seggeo->GetLayerCellGeom(layer);
    assert(layerGeom);

    if (nZBins <= 0)
    {
      nZBins = layerGeom->get_zbins();
      assert(nZBins > 0);
    }
    else
    {
      if ((int) nZBins != layerGeom->get_zbins())
      {
        cout << "TPCDataStreamEmulator::InitRun - Fatal Error - nZBin at layer " << layer << " is " << layerGeom->get_zbins()
             << ", which is different from previous layers of nZBin = " << nZBins << endl;
        exit(1);
      }
    }
  }  //   for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  if (Verbosity() >= VERBOSITY_SOME)
    cout << "TPCDataStreamEmulator::get_HistoManager - Making PHTFileServer " << m_outputFileNameBase + ".root"
         << endl;
  PHTFileServer::get().open(m_outputFileNameBase + ".root", "RECREATE");

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);

  TH1D* h = new TH1D("hNormalization",  //
                     "Normalization;Items;Summed quantity", 10, .5, 10.5);
  int i = 1;
  h->GetXaxis()->SetBinLabel(i++, "Event count");
  h->GetXaxis()->SetBinLabel(i++, "Collision count");
  h->GetXaxis()->SetBinLabel(i++, "TPC G4Hit");
  h->GetXaxis()->SetBinLabel(i++, "TPC G4Hit Edep");
  h->GetXaxis()->SetBinLabel(i++, "TPC Hit");
  h->GetXaxis()->SetBinLabel(i++, "TPC Wavelet");
  h->GetXaxis()->SetBinLabel(i++, "TPC DataSize");

  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  //  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  //  {
  //    const PHG4CylinderCellGeom* layer_geom = seggeo->GetLayerCellGeom(layer);

  //    const string histNameCellHit(boost::str(boost::format{"hCellHit_Layer%1%"} % layer));
  //    const string histNameCellCharge(boost::str(boost::format{"hCellCharge_Layer%1%"} % layer));

  //  }

  hm->registerHisto(m_hDataSize =
                        new TH1D("hDataSize",  //
                                 "TPC Data Size per Event;Data size [Byte];Count",
                                 10000, 0, 20e6));

  hm->registerHisto(m_hWavelet =
                        new TH1D("hWavelet",  //
                                 "TPC Recorded Wavelet per Event;Wavelet count;Count",
                                 10000, 0, 4e6));

  hm->registerHisto(m_hNChEta =
                        new TH1D("hNChEta",  //
                                 "Charged particle #eta distribution;#eta;Count",
                                 1000, -5, 5));

  hm->registerHisto(m_hLayerWaveletSize =
                        new TH2D("hLayerWaveletSize",  //
                                 "Number of Recorded ADC sample per Wavelet;Layer ID;ADC Sample Count per Wavelet",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 nZBins, -.5, nZBins - .5));

  hm->registerHisto(m_hLayerHit =
                        new TH2D("hLayerHit",  //
                                 "Number of Recorded ADC sample per channel;Layer ID;ADC Sample Count",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 nZBins, -.5, nZBins - .5));

  hm->registerHisto(m_hLayerDataSize =
                        new TH2D("hLayerDataSize",  //
                                 "Data size per channel;Layer ID;Data size [Byte]",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 2 * nZBins, 0, 2 * nZBins));

  hm->registerHisto(m_hLayerSumHit =
                        new TH2D("hLayerSumHit",  //
                                 "Number of Recorded ADC sample per layer;Layer ID;ADC Sample Count",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 10000, -.5, 99999.5));

  hm->registerHisto(m_hLayerSumDataSize =
                        new TH2D("hLayerSumDataSize",  //
                                 "Data size per trigger per layer;Layer ID;Data size [Byte]",
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
                                 1000, 0, .5e6));

  hm->registerHisto(m_hLayerZBinHit =
                        new TH2D("hLayerZBinHit",  //
                                 "Number of Recorded ADC sample per Time Bin;z bin ID;Layer ID",
                                 nZBins, -.5, nZBins - .5,
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5));

  hm->registerHisto(m_hLayerZBinADC =
                        new TH2D("hLayerZBinADC",  //
                                 "Sum ADC per Time Bin;z bin ID;Layer ID",
                                 nZBins, -.5, nZBins - .5,
                                 m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5));
  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCDataStreamEmulator::process_event(PHCompositeNode* topNode)
{
  m_evtCounter += 1;

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);
  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto("hNormalization"));
  assert(h_norm);
  h_norm->Fill("Event count", 1);

  PHG4HitContainer* g4hit = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_SVTX");
  if (!g4hit)
  {
    cout << "TPCDataStreamEmulator::process_event - Could not locate g4 hit node G4HIT_SVTX" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  SvtxHitMap* hits = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  if (!hits)
  {
    cout << "PCDataStreamEmulator::process_event - ERROR: Can't find node SvtxHitMap" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_SVTX");
  if (!cells)
  {
    cout << "TPCDataStreamEmulator::process_event - could not locate cell node "
         << "G4CELL_SVTX" << endl;
    exit(1);
  }

  PHG4CylinderCellGeomContainer* seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!seggeo)
  {
    cout << "TPCDataStreamEmulator::process_event - could not locate geo node "
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

      cout << "TPCDataStreamEmulator::process_event - - missing node PHHepMCGenEventMap. Skipping HepMC stat." << std::endl;
    }
  }
  else
  {
    h_norm->Fill("Collision count", geneventmap->size());
  }

  PHG4TruthInfoContainer* truthInfoList = findNode::getClass<PHG4TruthInfoContainer>(topNode,
                                                                                     "G4TruthInfo");
  if (!truthInfoList)
  {
    cout << "TPCDataStreamEmulator::process_event - Fatal Error - "
         << "unable to find DST node "
         << "G4TruthInfo" << endl;
    assert(truthInfoList);
  }

  PHG4TruthInfoContainer::ConstRange primary_range =
      truthInfoList->GetPrimaryParticleRange();

  for (PHG4TruthInfoContainer::ConstIterator particle_iter =
           primary_range.first;
       particle_iter != primary_range.second;
       ++particle_iter)
  {
    const PHG4Particle* p = particle_iter->second;
    assert(p);

    TParticlePDG* pdg_p = TDatabasePDG::Instance()->GetParticle(
        p->get_pid());
    assert(pdg_p);

    if (fabs(pdg_p->Charge()) > 0)
    {
      TVector3 pvec(p->get_px(), p->get_py(), p->get_pz());

      if (pvec.Perp2()>0)
      {
        assert(m_hNChEta);
        m_hNChEta->Fill(pvec.PseudoRapidity());
      }
    }

  }  //          if (_load_all_particle) else

  for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(layer);

    for (PHG4HitContainer::ConstIterator hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      const double edep = hiter->second->get_edep();
      h_norm->Fill("TPC G4Hit Edep", edep);
      h_norm->Fill("TPC G4Hit", 1);
    }  //     for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)

  }  //   for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  // prepreare stat. storage
  int nZBins = 0;
  vector<array<vector<int>, 2> > layerChanHit(m_maxLayer + 1);
  vector<array<vector<int>, 2> > layerChanDataSize(m_maxLayer + 1);
  int nWavelet = 0;
  int sumDataSize = 0;
  for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    PHG4CylinderCellGeom* layerGeom =
        seggeo->GetLayerCellGeom(layer);
    assert(layerGeom);

    // start with an empty vector of cells for each phibin
    const int nphibins = layerGeom->get_phibins();
    assert(nphibins > 0);

    if (Verbosity() >= VERBOSITY_MORE)
    {
      cout << "TPCDataStreamEmulator::process_event - init layer " << layer << " with "
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
        cout << "TPCDataStreamEmulator::process_event - Fatal Error - nZBin at layer " << layer << " is " << layerGeom->get_zbins()
             << ", which is different from previous layers of nZBin = " << nZBins << endl;
        exit(1);
      }
    }

    for (unsigned int side = 0; side < 2; ++side)
    {
      layerChanHit[layer][side].resize(nphibins, 0);

      layerChanDataSize[layer][side].resize(nphibins, 0);
    }  //     for (unsigned int side = 0; side < 2; ++side)

  }  //   for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  assert(nZBins > 0);

  // count hits and make wavelets
  int last_layer = -1;
  int last_side = -1;
  int last_phibin = -1;
  int last_zbin = -1;
  vector<unsigned int> last_wavelet;
  int last_wavelet_hittime = -1;

  for (SvtxHitMap::Iter iter = hits->begin(); iter != hits->end(); ++iter)
  {
    SvtxHit* hit = iter->second;

    const int layer = hit->get_layer();

    if (layer < m_minLayer or layer > m_maxLayer) continue;

    PHG4Cell* cell = cells->findCell(hit->get_cellid());                           //not needed once geofixed
    const int phibin = PHG4CellDefs::SizeBinning::get_phibin(cell->get_cellid());  //cell->get_binphi();
    const int zbin = PHG4CellDefs::SizeBinning::get_zbin(cell->get_cellid());      //cell->get_binz();
    const int side = (zbin < nZBins / 2) ? 0 : 1;

    // new wavelet?
    if (last_layer != layer or last_phibin != phibin or last_side != side or abs(last_zbin - zbin) != 1)
    {
      // save last wavelet
      if (last_wavelet.size() > 0)
      {
        const int datasize = writeWavelet(last_layer, last_side, last_phibin, last_wavelet_hittime, last_wavelet);
        assert(datasize > 0);

        nWavelet += 1;
        sumDataSize += datasize;
        layerChanDataSize[last_layer][last_side][last_phibin] += datasize;

        last_wavelet.clear();
        last_zbin = -1;
      }

      // z-R cut on digitized wavelet
      PHG4CylinderCellGeom* layerGeom =
          seggeo->GetLayerCellGeom(layer);
      assert(layerGeom);
      const double z_abs = fabs(layerGeom->get_zcenter(zbin));
      const double r = layerGeom->get_radius();
      TVector3 acceptanceVec(r, 0, z_abs - m_vertexZAcceptanceCut);
      const double eta = acceptanceVec.PseudoRapidity();

      if (eta > m_etaAcceptanceCut) continue;

      // make new wavelet
      last_layer = layer;
      last_side = side;
      last_phibin = phibin;

      // time check
      last_wavelet_hittime = (side == 0) ? (zbin) : (nZBins - 1 - zbin);
      assert(last_wavelet_hittime >= 0);
      assert(last_wavelet_hittime <= nZBins / 2);
    }  //     if (last_layer != layer or last_phibin != phibin)

    if (Verbosity() >= VERBOSITY_A_LOT)
    {
      cout << "TPCDataStreamEmulator::process_event -  layer " << layer << " hit with "

           << "phibin = " << phibin
           << ",zbin = " << zbin
           << ",side = " << side
           << ",last_wavelet.size() = " << last_wavelet.size()
           << ",last_zbin = " << last_zbin
           << endl;
    }

    // more checks on signal continuity
    if (last_wavelet.size() > 0)
    {
      if (side == 0)
      {
        assert(zbin - last_zbin == 1);
      }
      else
      {
        assert(last_zbin - zbin == 1);
      }
    }

    // record adc
    unsigned int adc = hit->get_adc();
    last_wavelet.push_back(adc);
    last_zbin = zbin;

    // statistics
    layerChanHit[layer][side][phibin] += 1;
    assert(m_hLayerZBinHit);
    m_hLayerZBinHit->Fill(zbin, layer, 1);
    assert(m_hLayerZBinADC);
    m_hLayerZBinADC->Fill(zbin, layer, adc);

  }  //   for(SvtxHitMap::Iter iter = hits->begin(); iter != hits->end(); ++iter) {

  // save last wavelet
  if (last_wavelet.size() > 0)
  {
    const int datasize = writeWavelet(last_layer, last_side, last_phibin, last_wavelet_hittime, last_wavelet);
    assert(datasize > 0);

    nWavelet += 1;
    sumDataSize += datasize;
    layerChanDataSize[last_layer][last_side][last_phibin] += datasize;
  }

  // statistics
  for (int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  {
    for (unsigned int side = 0; side < 2; ++side)
    {
      int sumHit = 0;
      for (const int& hit : layerChanHit[layer][side])
      {
        sumHit += hit;

        assert(m_hLayerHit);
        m_hLayerHit->Fill(layer, hit);
        h_norm->Fill("TPC Hit", hit);

        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << "TPCDataStreamEmulator::process_event - hit: "
               << "hit = " << hit
               << "at layer = " << layer
               << ", side = " << side
               << endl;
        }
      }

      if (Verbosity() >= VERBOSITY_MORE)
      {
        cout << "TPCDataStreamEmulator::process_event - hLayerSumCellHit->Fill(" << layer << ", " << sumHit << ")" << endl;
      }
      assert(m_hLayerSumHit);
      m_hLayerSumHit->Fill(layer, sumHit);

      double sumData = 0;
      for (const int& data : layerChanDataSize[layer][side])
      {
        sumData += data;

        assert(m_hLayerDataSize);
        m_hLayerDataSize->Fill(layer, data);
      }
      assert(m_hLayerSumDataSize);
      m_hLayerSumDataSize->Fill(layer, sumData);
    }  //    for (unsigned int side = 0; side < 2; ++side)

  }  //  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)

  assert(m_hWavelet);
  m_hWavelet->Fill(nWavelet);
  h_norm->Fill("TPC Wavelet", nWavelet);
  assert(m_hDataSize);
  m_hDataSize->Fill(sumDataSize);
  h_norm->Fill("TPC DataSize", sumDataSize);

  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCDataStreamEmulator::writeWavelet(int layer, int side, int phibin, int hittime, const vector<unsigned int>& wavelet)
{
  static const int headersize = 2;  // 2-byte header per wavelet

  //data in byte aligned and padding
  const int datasizebit = wavelet.size() * 10;
  int datasizebyte = datasizebit / 8;
  if (datasizebyte * 8 < datasizebit) datasizebyte += 1;

  assert(m_hLayerWaveletSize);
  m_hLayerWaveletSize->Fill(layer, wavelet.size());

  return headersize + datasizebyte;
}

Fun4AllHistoManager*
TPCDataStreamEmulator::getHistoManager()
{
  static string histname("TPCDataStreamEmulator_HISTOS");

  Fun4AllServer* se = Fun4AllServer::instance();
  Fun4AllHistoManager* hm = se->getHistoManager(histname);

  if (not hm)
  {
    cout
        << "TPCDataStreamEmulator::get_HistoManager - Making Fun4AllHistoManager " << histname
        << endl;
    hm = new Fun4AllHistoManager(histname);
    se->registerHistoManager(hm);
  }

  assert(hm);

  return hm;
}
