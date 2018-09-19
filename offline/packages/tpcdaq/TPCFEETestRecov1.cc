/*
 * TPCFEETestRecov1.cc
 *
 *  Created on: Sep 19, 2018
 *      Author: jinhuang
 */

#include "TPCFEETestRecov1.h"

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

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>
#include <Event/packetConstants.h>

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
using namespace TPCDaqDefs::FEEv1;

TPCFEETestRecov1::TPCFEETestRecov1(const std::string& outputfilename)
  : SubsysReco("TPCFEETestRecov1")
  , m_outputFileName(outputfilename)
  , m_T(nullptr)
  , m_event(-1)
{
}

TPCFEETestRecov1::~TPCFEETestRecov1()
{
}

int TPCFEETestRecov1::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCFEETestRecov1::End(PHCompositeNode* topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
    cout << "TPCFEETestRecov1::End - write to " << m_outputFileName << endl;
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

int TPCFEETestRecov1::InitRun(PHCompositeNode* topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
    cout << "TPCFEETestRecov1::get_HistoManager - Making PHTFileServer " << m_outputFileName
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

  m_T = new TTree("T", "TPC FEE Tree");
  assert(m_T);

  //  for (unsigned int layer = m_minLayer; layer <= m_maxLayer; ++layer)
  //  {
  //    const PHG4CylinderCellGeom* layer_geom = seggeo->GetLayerCellGeom(layer);

  //    const string histNameCellHit(boost::str(boost::format{"hCellHit_Layer%1%"} % layer));
  //    const string histNameCellCharge(boost::str(boost::format{"hCellCharge_Layer%1%"} % layer));

  //  }

  //  hm->registerHisto(new TH2D("hLayerCellHit",  //
  //                             "Number of ADC time-bin hit per channel;Layer ID;Hit number",
  //                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
  //                             300, -.5, 299.5));
  //  hm->registerHisto(new TH2D("hLayerCellCharge",  //
  //                             "Charge integrated over drift window per channel;Layer ID;Charge [fC]",
  //                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
  //                             1000, 0, 1e7 * eplus / (1e-15 * coulomb)));
  //
  //  hm->registerHisto(new TH2D("hLayerSumCellHit",  //
  //                             "Number of ADC time-bin hit integrated over channels per layer;Layer ID;Hit number",
  //                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
  //                             10000, -.5, 99999.5));
  //  hm->registerHisto(new TH2D("hLayerSumCellCharge",  //
  //                             "Charge integrated over drift window and channel per layer;Layer ID;Charge [fC]",
  //                             m_maxLayer - m_minLayer + 1, m_minLayer - .5, m_maxLayer + .5,
  //                             10000, 0, 1000 * 4e6 * eplus / (1e-15 * coulomb)));

  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCFEETestRecov1::process_event(PHCompositeNode* topNode)
{
  ++m_event;

  Fun4AllHistoManager* hm = getHistoManager();
  assert(hm);
  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto("hNormalization"));
  assert(h_norm);

  Event* event = findNode::getClass<Event>(topNode, "PRDF");
  if (event == nullptr)
  {
    if (Verbosity() >= VERBOSITY_SOME)
      cout << "GenericUnpackPRDF::Process_Event - Event not found" << endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  if (verbosity >= VERBOSITY_SOME)
    event->identify();

  // search for data event
  if (event->getEvtType() != DATAEVENT)
    return Fun4AllReturnCodes::DISCARDEVENT;

  Packet* p = event->getPacket(1024, ID4EVT);
  if (p == nullptr)
    return Fun4AllReturnCodes::DISCARDEVENT;

  if (verbosity >= VERBOSITY_MORE)
    p->identify();

  uint32_t bx_seed = 0;
  for (int channel = 0; channel < N_CHANNELS; channel++)
  {
    uint16_t size = p->iValue(channel * kPACKET_LENGTH + 1) & 0xffff;        // number of words until the next channel (header included). this is the real packet_length
    uint8_t packet_type = p->iValue(channel * kPACKET_LENGTH + 2) & 0xffff;  // that's the Elink packet type
    uint32_t bx_counter = ((p->iValue(channel * kPACKET_LENGTH + 4) & 0xffff) << 4) | (p->iValue(channel * kPACKET_LENGTH + 5) & 0xffff);
    uint8_t sampa_address = (p->iValue(channel * kPACKET_LENGTH + 3) >> 5) & 0xf;
    uint16_t sampa_channel = p->iValue(channel * kPACKET_LENGTH + 3) & 0x1f;
    uint16_t fee_channel = (sampa_address << 5) | sampa_channel;

    if (channel == 0)
    {
      bx_seed = bx_counter;
    }
    else if (bx_seed != bx_counter)
    {
      printf("TPCFEETestRecov1::process_event - ERROR: Malformed packet, event number %i, reason: bx_counter mismatch (expected 0x%x, got 0x%x)\n", m_event, bx_seed, bx_counter);

      event->identify();
      p->identify();
      return Fun4AllReturnCodes::DISCARDEVENT;
    }

    if (fee_channel > 255 || sampa_address > 7 || sampa_channel > 31)
    {
      printf("TPCFEETestRecov1::process_event - ERROR: Malformed packet, event number %i, reason: bad channel (got %i, sampa_addr: %i, sampa_chan: %i)\n", m_event, fee_channel, sampa_address, sampa_channel);

      event->identify();
      p->identify();
      return Fun4AllReturnCodes::DISCARDEVENT;
    }

    //    SampaChannel *chan = fee_data->append(new SampaChannel(fee_channel, bx_counter, packet_type));

    for (int sample = 0; sample < kSAMPLE_LENGTH; sample++)
    {
      //        chan->append(p->iValue(channel * PACKET_LENGTH + 9 + sample) & 0xffff);
    }
  }

  h_norm->Fill("Event count", 1);
  m_T->Fill();

  return Fun4AllReturnCodes::EVENT_OK;
}
Fun4AllHistoManager*
TPCFEETestRecov1::getHistoManager()
{
  static string histname("TPCFEETestRecov1_HISTOS");

  Fun4AllServer* se = Fun4AllServer::instance();
  Fun4AllHistoManager* hm = se->getHistoManager(histname);

  if (not hm)
  {
    cout
        << "TPCFEETestRecov1::get_HistoManager - Making Fun4AllHistoManager "
        << histname << endl;
    hm = new Fun4AllHistoManager(histname);
    se->registerHistoManager(hm);
  }

  assert(hm);

  return hm;
}
