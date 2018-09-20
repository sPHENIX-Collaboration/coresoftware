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
//#include <Event/packetConstants.h>
#include <Event/oncsSubConstants.h>

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
  , m_eventT(nullptr)
  , m_chanT(nullptr)
  , m_event(-1)
  , m_pchanHeader(&m_chanHeader)
  , m_chanData(kSAMPLE_LENGTH, 0)
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

  m_eventT->Write();
  m_chanT->Write();

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

  m_eventT = new TTree("eventT", "TPC FEE per-event Tree");
  assert(m_eventT);
  m_eventT->Branch("event", &m_event, "event/I");

  m_chanT = new TTree("chanT", "TPC FEE per-channel Tree");
  assert(m_chanT);
  m_chanT->Branch("event", &m_event, "event/I");
  m_chanT->Branch("header", &m_pchanHeader);
  m_chanT->Branch("adc", m_chanData.data(), str(boost::format("adc[%d]/i") % kSAMPLE_LENGTH).c_str());

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

  m_event = event->getEvtSequence();

  Packet* p = event->getPacket(kPACKET_ID, ID4EVT);
  if (p == nullptr)
    return Fun4AllReturnCodes::DISCARDEVENT;

  if (verbosity >= VERBOSITY_SOME) p->identify();

  if (verbosity >= VERBOSITY_MORE)
  {
    cout << "TPCFEETestRecov1::process_event - p->iValue(0) = "
         << p->iValue(0) << ", p->iValue(1) = " << p->iValue(1)
         << ", p->iValue(2) = " << p->iValue(2)
         << ", p->iValue(3) = " << p->iValue(3) << endl;
    p->dump();
  }

  uint32_t bx_seed = 0;
  for (unsigned int channel = 0; channel < kN_CHANNELS; channel++)
  {
    m_chanHeader.m_size = p->iValue(channel * kPACKET_LENGTH + 1) & 0xffff;         // number of words until the next channel (header included). this is the real packet_length
    m_chanHeader.m_packet_type = p->iValue(channel * kPACKET_LENGTH + 2) & 0xffff;  // that's the Elink packet type
    m_chanHeader.m_bx_counter = ((p->iValue(channel * kPACKET_LENGTH + 4) & 0xffff) << 4) | (p->iValue(channel * kPACKET_LENGTH + 5) & 0xffff);
    m_chanHeader.m_sampa_address = (p->iValue(channel * kPACKET_LENGTH + 3) >> 5) & 0xf;
    m_chanHeader.m_sampa_channel = p->iValue(channel * kPACKET_LENGTH + 3) & 0x1f;
    m_chanHeader.m_fee_channel = (m_chanHeader.m_sampa_address << 5) | m_chanHeader.m_sampa_channel;

    const pair<int, int> pad = SAMPAChan2PadXY(m_chanHeader.m_fee_channel);

    m_chanHeader.m_pad_x = pad.first;
    m_chanHeader.m_pad_y = pad.second;

    if (channel == 0)
    {
      bx_seed = m_chanHeader.m_bx_counter;
    }
    else if (bx_seed != m_chanHeader.m_bx_counter)
    {
      printf("TPCFEETestRecov1::process_event - ERROR: Malformed packet, event number %i, reason: bx_counter mismatch (expected 0x%x, got 0x%x)\n", m_event, bx_seed, m_chanHeader.m_bx_counter);

      event->identify();
      p->identify();
      return Fun4AllReturnCodes::DISCARDEVENT;
    }

    if (m_chanHeader.m_fee_channel > 255 || m_chanHeader.m_sampa_address > 7 || m_chanHeader.m_sampa_channel > 31)
    {
      printf("TPCFEETestRecov1::process_event - ERROR: Malformed packet, event number %i, reason: bad channel (got %i, sampa_addr: %i, sampa_chan: %i)\n", m_event, m_chanHeader.m_fee_channel, m_chanHeader.m_sampa_address, m_chanHeader.m_sampa_channel);

      event->identify();
      p->identify();
      return Fun4AllReturnCodes::DISCARDEVENT;
    }

    //    SampaChannel *chan = fee_data->append(new SampaChannel(fee_channel, bx_counter, packet_type));

    assert(m_chanData.size() == kSAMPLE_LENGTH);
    fill(m_chanData.begin(), m_chanData.end(), 0);
    for (unsigned int sample = 0; sample < kSAMPLE_LENGTH; sample++)
    {
      //        chan->append(p->iValue(channel * PACKET_LENGTH + 9 + sample) & 0xffff);
      uint32_t value = p->iValue(channel * kPACKET_LENGTH + 9 + sample) & 0xffff;
      m_chanData[sample] = value;
    }

    if (verbosity >= VERBOSITY_MORE)
    {
      cout << "TPCFEETestRecov1::process_event - "
           << "m_chanHeader.m_size = " << int(m_chanHeader.m_size) << ", "
           << "m_chanHeader.m_packet_type = " << int(m_chanHeader.m_packet_type) << ", "
           << "m_chanHeader.m_bx_counter = " << int(m_chanHeader.m_bx_counter) << ", "
           << "m_chanHeader.m_sampa_address = " << int(m_chanHeader.m_sampa_address) << ", "
           << "m_chanHeader.m_sampa_channel = " << int(m_chanHeader.m_sampa_channel) << ", "
           << "m_chanHeader.m_fee_channel = " << int(m_chanHeader.m_fee_channel) << ": "
           << " ";

      for (unsigned int sample = 0; sample < kSAMPLE_LENGTH; sample++)
      {
        cout << "data[" << sample << "] = " << m_chanData[sample] << " ";
      }

      cout << endl;
    }

    m_chanT->Fill();
  }

  h_norm->Fill("Event count", 1);
  m_eventT->Fill();

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
