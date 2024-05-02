// Include necessary files
#include "TpcChanQA.h"

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
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <string>
//

//____________________________________________________________________________..
TpcChanQA::TpcChanQA(const std::string &name)
  : SubsysReco("TpcChanQA")
  , m_fname(name)
{
  // reserves memory for max ADC samples
  m_adcSamples.resize(1024, 0);
}

//____________________________________________________________________________..
int TpcChanQA::InitRun(PHCompositeNode * /*unused*/)
{
  // Takes string of raw data file and truncates it down to sector number
  sectorNum = m_fname;
  size_t pos = sectorNum.find("TPC_ebdc");
  sectorNum.erase(sectorNum.begin(), sectorNum.begin() + pos + 8);
  sectorNum.erase(sectorNum.begin() + 2, sectorNum.end());
  if (sectorNum.at(0) == '0')
  {
    sectorNum.erase(sectorNum.begin(), sectorNum.begin() + 1);
  }
  //

  // Sets side to South if SectorNum > 11 (EBDC 12-23)
  if (stoi(sectorNum) > 11)
  {
    side = 1;
  }

  // Creates data file and checks whether it was successfully opened
  m_file = TFile::Open(m_fname.c_str(), "recreate");
  assert(m_file->IsOpen());

  // Define histograms initialized in header file
  std::string name = "h_channel_hits_sec" + sectorNum;
  h_channel_hits = new TH1F(name.c_str(), name.c_str(), 256, 0, 256);
  name = "h_channel_ADCs_sec" + sectorNum;
  h_channel_ADCs = new TH2F(name.c_str(), name.c_str(), 256, 0, 256, 1024, 0, 1024);
  //

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcChanQA::process_event(PHCompositeNode *topNode)
{
  // Defines object from class Event which calls getClass function from
  // findNode class
  Event *_event = findNode::getClass<Event>(topNode, "PRDF");

  // Checks if event exists and returns error if not
  if (_event == nullptr)
  {
    std::cout << "TPCRawDataTree::Process_Event - Event not found" << std::endl;
    return -1;
  }

  // Checks if event is "special" and discards it if so
  if (_event->getEvtType() >= 8)  /// special events
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  // Loop over packets in event
  for (int packet : m_packets)
  {
    if (Verbosity())
    {
      std::cout << __PRETTY_FUNCTION__ << " : decoding packet " << packet << std::endl;
    }

    m_packet = packet;

    // Assigns data from packet to variable p then checks if packet exists
    // Continues if not
    std::unique_ptr<Packet> p(_event->getPacket(m_packet));
    if (!p)
    {
      if (Verbosity())
      {
        std::cout << __PRETTY_FUNCTION__ << " : missing packet " << packet << std::endl;
      }
      continue;
    }

    // pull number of waveforms
    m_nWaveformInFrame = p->iValue(0, "NR_WF");

    for (int wf = 0; wf < m_nWaveformInFrame; wf++)
    {
      m_Channel = p->iValue(wf, "CHANNEL");
      m_nSamples = p->iValue(wf, "SAMPLES");

      h_channel_hits->Fill(m_Channel);

      // Checks if sample number and number of ADC values agrees
      assert(m_nSamples < (int) m_adcSamples.size());

      // Loop over samples in waveform
      for (int s = 0; s < m_nSamples; s++)
      {
        // Assign ADC value of sample s in waveform wf to adcSamples[s]
        m_adcSamples[s] = p->iValue(wf, s);
        h_channel_ADCs->Fill(m_Channel, m_adcSamples[s]);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcChanQA::End(PHCompositeNode * /*unused*/)
{
  // Set histogram directory to 0 so data is saved after closing file
  h_channel_hits->SetDirectory(nullptr);
  h_channel_ADCs->SetDirectory(nullptr);

  // Write histograms to file
  m_file->cd();
  h_channel_hits->Write();
  h_channel_ADCs->Write();

  std::cout << __PRETTY_FUNCTION__ << " : completed saving to " << m_file->GetName() << std::endl;
  m_file->ls();

  // Close the file
  m_file->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}
