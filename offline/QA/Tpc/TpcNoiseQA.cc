// Includes
#include "TpcNoiseQA.h"

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/packet.h>

#include <TH1.h>
#include <TH2.h>

#include <boost/format.hpp>

#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
//

//____________________________________________________________________________..
TpcNoiseQA::TpcNoiseQA(const std::string &name)
  : SubsysReco("TpcNoiseQA")
  , m_fname(name)
{
  // reserves memory for max ADC samples
  m_adcSamples.resize(1024, 0);

  M.setMapNames("AutoPad-R1-RevA.sch.ChannelMapping.csv", "AutoPad-R2-RevA-Pads.sch.ChannelMapping.csv", "AutoPad-R3-RevA.sch.ChannelMapping.csv");
}

//____________________________________________________________________________..
int TpcNoiseQA::InitRun(PHCompositeNode * /*unused*/)
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

  // double r_bins_new[r_bins_N + 1];
  for (int i = 0; i < r_bins_N + 1; i++)
  {
    r_bins_new[i] = -1 * r_bins[r_bins_N - i];
    r_bins_new[i + r_bins_N + nphi + 1] = r_bins[i];
  }
  for (int i = 0; i < nphi + 1; i++)
  {
    r_bins_new[i + r_bins_N] = phi_bins[i];
  }

  createHistos();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcNoiseQA::process_event(PHCompositeNode *topNode)
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

  for (int fee_no = 0; fee_no < 26; fee_no++)
  {
    for (int channel_no = 0; channel_no < 256; channel_no++)
    {
      ave_adc_fee_channel[fee_no][channel_no] = 0.0;
      std_adc_fee_channel[fee_no][channel_no] = 0.0;
      counts_adc_fee_channel[fee_no][channel_no] = 0.0;
    }
  }

  // Call HistoManager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // Reference histograms initialized in header file to histos in HistoManager
  h_NPol_Ped_Mean = dynamic_cast<TH2F *>(hm->getHisto(boost::str(boost::format("%sNPol_Ped_Mean_sec%s") % getHistoPrefix() % sectorNum.c_str()).c_str()));
  h_NPol_Ped_RMS = dynamic_cast<TH2F *>(hm->getHisto(boost::str(boost::format("%sNPol_Ped_RMS_sec%s") % getHistoPrefix() % sectorNum.c_str()).c_str()));
  h_SPol_Ped_Mean = dynamic_cast<TH2F *>(hm->getHisto(boost::str(boost::format("%sSPol_Ped_Mean_sec%s") % getHistoPrefix() % sectorNum.c_str()).c_str()));
  h_SPol_Ped_RMS = dynamic_cast<TH2F *>(hm->getHisto(boost::str(boost::format("%sSPol_Ped_RMS_sec%s") % getHistoPrefix() % sectorNum.c_str()).c_str()));
  //

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
      m_FEE = p->iValue(wf, "FEE");
      m_Channel = p->iValue(wf, "CHANNEL");
      m_nSamples = p->iValue(wf, "SAMPLES");

      // Checks if sample number and number of ADC values agrees
      // assert(m_nSamples < (int) m_adcSamples.size());
      if (m_nSamples > (int) m_adcSamples.size())
      {
        continue;
      }

      dead = false;
      // Loop over samples in waveform
      for (int s = 0; s < m_nSamples; s++)
      {
        // Assign ADC value of sample s in waveform wf to adcSamples[s]
        m_adcSamples[s] = p->iValue(wf, s);

        if (m_adcSamples[s] == 0 || std::isnan(float(m_adcSamples[s])))
        {
          dead = true;
          break;
        }
      }

      if (dead)
      {
        continue;
      }

      for (int adc_sam_no = 0; adc_sam_no < m_nSamples; adc_sam_no++)
      {
        if (m_adcSamples[adc_sam_no] < 1024)
        {
          ave_adc_fee_channel[m_FEE][m_Channel] += m_adcSamples[adc_sam_no];
          std_adc_fee_channel[m_FEE][m_Channel] += pow(m_adcSamples[adc_sam_no], 2);
          counts_adc_fee_channel[m_FEE][m_Channel] += 1.0;
        }
      }
    }
  }

  for (int fee_no = 0; fee_no < 26; fee_no++)
  {
    for (int channel_no = 0; channel_no < 256; channel_no++)
    {
      if (counts_adc_fee_channel[fee_no][channel_no] != 0.0)
      {
        temp1 = ave_adc_fee_channel[fee_no][channel_no] / counts_adc_fee_channel[fee_no][channel_no];
        temp2 = std_adc_fee_channel[fee_no][channel_no] / counts_adc_fee_channel[fee_no][channel_no];
        ave_adc_fee_channel[fee_no][channel_no] = temp1;
        std_adc_fee_channel[fee_no][channel_no] = temp2;
      }
      else
      {
        ave_adc_fee_channel[fee_no][channel_no] = 0.0;
        std_adc_fee_channel[fee_no][channel_no] = 0.0;
      }

      // setting the mapp of the FEE
      feeM = FEE_map[fee_no];
      if (mod_arr[fee_no] == 2)
      {
        feeM += 6;
      }
      if (mod_arr[fee_no] == 3)
      {
        feeM += 14;
      }
      key = 256 * (feeM) + channel_no;
      R = M.getR(feeM, channel_no);
      phi = pow(-1, side) * M.getPhi(feeM, channel_no) + (stoi(sectorNum) - side * 12.0) * M_PI / 6.0;

      if (phi < 0.0)
      {
        phi = phi + 2.0 * M_PI;
      }

      pedMean = ave_adc_fee_channel[fee_no][channel_no];
      pedStdi = sqrt(std_adc_fee_channel[fee_no][channel_no] - pow(ave_adc_fee_channel[fee_no][channel_no], 2));

      if (side == 0)
      {
        h_NPol_Ped_Mean->Fill(phi, R, pedMean);
        h_NPol_Ped_RMS->Fill(phi, R, pedStdi);
      }

      else
      {
        h_SPol_Ped_Mean->Fill(phi, R, pedMean);
        h_SPol_Ped_RMS->Fill(phi, R, pedStdi);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcNoiseQA::End(PHCompositeNode * /*unused*/)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
std::string TpcNoiseQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }  // Define prefix to all histos in HistoManager

//____________________________________________________________________________..
void TpcNoiseQA::createHistos()
{
  // Initialize HistoManager
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // Create and register histos in HistoManager
  {
    auto h = new TH2F(boost::str(boost::format("%sNPol_Ped_Mean_sec%s") % getHistoPrefix() % sectorNum.c_str()).c_str(), ";x;y", (2 * r_bins_N + nphi + 1), r_bins_new, (2 * r_bins_N + nphi + 1), r_bins_new);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(boost::str(boost::format("%sNPol_Ped_RMS_sec%s") % getHistoPrefix() % sectorNum.c_str()).c_str(), ";x;y", (2 * r_bins_N + nphi + 1), r_bins_new, (2 * r_bins_N + nphi + 1), r_bins_new);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(boost::str(boost::format("%sSPol_Ped_Mean_sec%s") % getHistoPrefix() % sectorNum.c_str()).c_str(), ";x;y", (2 * r_bins_N + nphi + 1), r_bins_new, (2 * r_bins_N + nphi + 1), r_bins_new);
    hm->registerHisto(h);
  }

  {
    auto h = new TH2F(boost::str(boost::format("%sSPol_Ped_RMS_sec%s") % getHistoPrefix() % sectorNum.c_str()).c_str(), ";x;y", (2 * r_bins_N + nphi + 1), r_bins_new, (2 * r_bins_N + nphi + 1), r_bins_new);
    hm->registerHisto(h);
  }
}
