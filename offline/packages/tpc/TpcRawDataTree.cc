
#include "TpcRawDataTree.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TFile.h>
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

  // add all possible TPC packet that we need to analyze
  // if a subset is needed,  call RemoveAllPackets()
  for (int packet = 4001; packet <= 4231; packet += 10)
  {
    m_packets.push_back(packet);
    m_packets.push_back(packet+1);
  }
}

//____________________________________________________________________________..
int TpcRawDataTree::InitRun(PHCompositeNode *)
{
  m_file = TFile::Open(m_fname.c_str(), "recreate");
  assert(m_file->IsOpen());

  m_PacketTree = new TTree("PacketTree", "Each entry is one packet");

  m_PacketTree->Branch("packet", &m_packet, "packet/I");
  m_PacketTree->Branch("frame", &m_frame, "frame/I");
  m_PacketTree->Branch("nWaveormInFrame", &m_nWaveormInFrame, "nWaveormInFrame/I");
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
    m_maxFEECount = p->iValue(0, "MAX_FEECOUNT");

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

      assert(m_nSamples < (int) m_adcSamples.size());  // no need for movements in memory allocation
      for (int s = 0; s < m_nSamples; s++)
      {
        m_adcSamples[s] = p->iValue(wf, s);
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
  m_file->Write();

  std::cout << __PRETTY_FUNCTION__ << " : completed saving to " << m_file->GetName() << std::endl;
  m_file->ls();

  m_file->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}
