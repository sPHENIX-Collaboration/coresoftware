
#include "TpcRawDataTree.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/packet.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include <cassert>
#include <iostream>
#include <memory>

//____________________________________________________________________________..
TpcRawDataTree::TpcRawDataTree(const std::string &name)
  : SubsysReco("TpcRawDataTree")
  , m_fname(name)
{
  // reserve memory for max ADC samples
  m_adcSamples.resize(1024, 0);
  M.setMapNames("AutoPad-R1-RevA.sch.ChannelMapping.csv", "AutoPad-R2-RevA-Pads.sch.ChannelMapping.csv", "AutoPad-R3-RevA.sch.ChannelMapping.csv");
}

//____________________________________________________________________________..
int TpcRawDataTree::InitRun(PHCompositeNode * /*unused*/)
{
  sectorNum = m_fname;
  size_t pos = sectorNum.find("TPC_ebdc");
  sectorNum.erase(sectorNum.begin(), sectorNum.begin() + pos + 8);
  sectorNum.erase(sectorNum.begin() + 2, sectorNum.end());
  if (sectorNum.at(0) == '0')
  {
    sectorNum.erase(sectorNum.begin(), sectorNum.begin() + 1);
  }
  if (stoi(sectorNum) > 11)
  {
    side = 1;
  }

  m_file = TFile::Open(m_fname.c_str(), "recreate");
  assert(m_file->IsOpen());

  m_PacketTree = new TTree("PacketTree", "Each entry is one packet");

  m_PacketTree->Branch("packet", &m_packet, "packet/I");
  m_PacketTree->Branch("frame", &m_frame, "frame/I");
  m_PacketTree->Branch("nWaveormInFrame", &m_nWaveormInFrame, "nWaveormInFrame/I");
  m_PacketTree->Branch("nTaggerInFrame", &m_nTaggerInFrame, "nTaggerInFrame/I");
  m_PacketTree->Branch("maxFEECount", &m_maxFEECount, "maxFEECount/I");

  m_SampleTree = new TTree("SampleTree", "Each entry is one waveform");

  m_SampleTree->Branch("packet", &m_packet, "packet/I");
  m_SampleTree->Branch("frame", &m_frame, "frame/I");
  m_SampleTree->Branch("nWaveormInFrame", &m_nWaveormInFrame, "nWaveormInFrame/I");
  m_SampleTree->Branch("maxFEECount", &m_maxFEECount, "maxFEECount/I");
  m_SampleTree->Branch("nSamples", &m_nSamples, "nSamples/I");
  m_SampleTree->Branch("adcSamples", &m_adcSamples[0], "adcSamples[nSamples]/s");
  m_SampleTree->Branch("fee", &m_fee, "fee/I");
  m_SampleTree->Branch("sampaAddress", &m_sampaAddress, "sampaAddress/I");
  m_SampleTree->Branch("sampaChannel", &m_sampaChannel, "sampaChannel/I");
  m_SampleTree->Branch("Channel", &m_Channel, "Channel/I");
  m_SampleTree->Branch("BCO", &m_BCO, "BCO/I");
  m_SampleTree->Branch("checksum", &m_checksum, "checksum/I");
  m_SampleTree->Branch("checksumError", &m_checksumError, "checksumError/I");

  m_TaggerTree = new TTree("TaggerTree", "Each entry is one tagger for level 1 trigger or endat tag");

  m_TaggerTree->Branch("packet", &m_packet, "packet/I");
  m_TaggerTree->Branch("frame", &m_frame, "frame/I");
  m_TaggerTree->Branch("tagger_type", &m_tagger_type, "tagger_type/s");
  m_TaggerTree->Branch("is_endat", &m_is_endat, "is_endat/b");
  m_TaggerTree->Branch("is_lvl1", &m_is_lvl1, "is_lvl1/b");
  m_TaggerTree->Branch("bco", &m_bco, "bco/l");
  m_TaggerTree->Branch("lvl1_count", &m_lvl1_count, "lvl1_count/i");
  m_TaggerTree->Branch("endat_count", &m_endat_count, "endat_count/i");
  m_TaggerTree->Branch("last_bco", &m_last_bco, "last_bco/l");
  m_TaggerTree->Branch("modebits", &m_modebits, "modebits/b");

  R1_hist = new TH1F("R1_hist", "R1_hist", 1024, -0.5, 1023.5);
  R2_hist = new TH1F("R2_hist", "R2_hist", 1024, -0.5, 1023.5);
  R3_hist = new TH1F("R3_hist", "R3_hist", 1024, -0.5, 1023.5);

  R1_time = new TH2F("R1_time", "R1_time", 360, -0.5, 359.5, 1024, -0.5, 1023.5);
  R2_time = new TH2F("R2_time", "R2_time", 360, -0.5, 359.5, 1024, -0.5, 1023.5);
  R3_time = new TH2F("R3_time", "R3_time", 360, -0.5, 359.5, 1024, -0.5, 1023.5);

  TotalFEE = new TH1F("TotalFEE", "Total FEE", 26, -0.5, 25.5);
  TotalFEEsampa = new TH1F("TotalFEEsampa", "Total FEE + sampa", 26 * 8, -0.5, 25 * 8 - .5);
  TotalFRAME = new TH1F("TotalFRAME", "Total FRAME", 21, -0.5, 20.5);

  checksumError_fee = new TH1F("FEEWithError", "FEE with Error", 26, -0.5, 25.5);
  checksumError_feesampa = new TH1F("FEEsampaWithError", "FEE*8+sampa with Error", 26 * 8, -0.5, 25 * 8 - .5);
  checksumError_frame = new TH1F("FRAMEWithError", "FRAME with Error", 21, -0.5, 20.5);

  if (m_includeXYPos)
  {
    m_SampleTree->Branch("xPos", &m_xPos, "xPos/d");
    m_SampleTree->Branch("yPos", &m_yPos, "yPos/d");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRawDataTree::process_event(PHCompositeNode *topNode)
{
  Event *_event = findNode::getClass<Event>(topNode, "PRDF");
  if (_event == nullptr)
  {
    std::cout << "TpcRawDataTree::Process_Event Event not found" << std::endl;
    return -1;
  }
  if (_event->getEvtType() >= 8)  /// special events
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  m_frame = _event->getEvtSequence();

  for (int packet : m_packets)
  {
    if (Verbosity())
    {
      std::cout << __PRETTY_FUNCTION__ << " : decoding packet " << packet << std::endl;
    }

    m_packet = packet;

    std::unique_ptr<Packet> p(_event->getPacket(m_packet));
    if (!p)
    {
      if (Verbosity())
      {
        std::cout << __PRETTY_FUNCTION__ << " : missing packet " << packet << std::endl;
      }

      continue;
    }

    m_nWaveormInFrame = p->iValue(0, "NR_WF");
    m_nTaggerInFrame = p->lValue(0, "N_TAGGER");
    m_maxFEECount = p->iValue(0, "MAX_FEECOUNT");

    for (int t = 0; t < m_nTaggerInFrame; t++)
    {
      /*uint16_t*/ m_tagger_type = (uint16_t) (p->lValue(t, "TAGGER_TYPE"));
      /*uint8_t*/ m_is_endat = (uint8_t) (p->lValue(t, "IS_ENDAT"));
      /*uint8_t*/ m_is_lvl1 = (uint8_t) (p->lValue(t, "IS_LEVEL1_TRIGGER"));
      /*uint64_t*/ m_bco = (uint64_t) (p->lValue(t, "BCO"));
      /*uint32_t*/ m_lvl1_count = (uint32_t) (p->lValue(t, "LEVEL1_COUNT"));
      /*uint32_t*/ m_endat_count = (uint32_t) (p->lValue(t, "ENDAT_COUNT"));
      /*uint64_t*/ m_last_bco = (uint64_t) (p->lValue(t, "LAST_BCO"));
      /*uint8_t*/ m_modebits = (uint8_t) (p->lValue(t, "MODEBITS"));

      m_TaggerTree->Fill();
    }

    for (int wf = 0; wf < m_nWaveormInFrame; wf++)
    {
      m_BCO = p->iValue(wf, "BCO");
      m_nSamples = p->iValue(wf, "SAMPLES");
      m_fee = p->iValue(wf, "FEE");

      m_sampaAddress = p->iValue(wf, "SAMPAADDRESS");
      m_sampaChannel = p->iValue(wf, "SAMPACHANNEL");
      m_Channel = p->iValue(wf, "CHANNEL");
      m_checksum = p->iValue(wf, "CHECKSUM");
      m_checksumError = p->iValue(wf, "CHECKSUMERROR");

      TH1 *fillHist;
      TH2 *fillHist2D;

      if (m_fee == 2 ||
          m_fee == 4 ||
          m_fee == 3 ||
          m_fee == 13 ||
          m_fee == 17 ||
          m_fee == 16)
      {
        fillHist = R1_hist;
        fillHist2D = R1_time;
      }
      else if (m_fee == 11 ||
               m_fee == 12 ||
               m_fee == 19 ||
               m_fee == 18 ||
               m_fee == 01 ||
               m_fee == 00 ||
               m_fee == 16 ||
               m_fee == 15)
      {
        fillHist = R2_hist;
        fillHist2D = R2_time;
      }
      else
      {
        fillHist = R3_hist;
        fillHist2D = R3_time;
      }

      assert(m_nSamples < (int) m_adcSamples.size());  // no need for movements in memory allocation
      for (int s = 0; s < m_nSamples; s++)
      {
        m_adcSamples[s] = p->iValue(wf, s);
        if (m_checksumError == 0)
        {
          fillHist->Fill(m_adcSamples[s]);
          fillHist2D->Fill(s, m_adcSamples[s]);
        }
        else
        {
          checksumError_fee->Fill(m_fee);
          checksumError_feesampa->Fill((m_fee * 8. + m_sampaAddress));
          checksumError_frame->Fill(m_frame);
        }
        TotalFEE->Fill(m_fee);
        TotalFEEsampa->Fill((m_fee * 8. + m_sampaAddress));
        TotalFRAME->Fill(m_frame);
      }
      if (m_includeXYPos)
      {
        int feeM = FEE_map[m_fee];
        if (FEE_R[m_fee] == 2)
        {
          feeM += 6;
        }
        if (FEE_R[m_fee] == 3)
        {
          feeM += 14;
        }
        int layer = M.getLayer(feeM, m_Channel);
        if (layer != 0)
        {
          double R = M.getR(feeM, m_Channel);
          double phi = M.getPhi(feeM, m_Channel) + (stod(sectorNum) - side * 12.) * M_PI / 6.;
          R /= 10.;  // convert mm to cm

          m_xPos = R * cos(phi);
          m_yPos = R * sin(phi);
        }
        else
        {
          m_xPos = 0.;
          m_yPos = 0.;
        }
      }
      m_SampleTree->Fill();
    }

    m_PacketTree->Fill();
  }  //   for (int packet : m_packets)

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcRawDataTree::End(PHCompositeNode * /*topNode*/)
{
  checksumError_fee->Divide(TotalFEE);
  checksumError_feesampa->Divide(TotalFEEsampa);
  checksumError_frame->Divide(TotalFRAME);

  TotalFEE->SetDirectory(nullptr);
  TotalFEEsampa->SetDirectory(nullptr);
  TotalFRAME->SetDirectory(nullptr);

  m_file->Write();

  std::cout << __PRETTY_FUNCTION__ << " : completed saving to " << m_file->GetName() << std::endl;
  m_file->ls();

  m_file->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}
