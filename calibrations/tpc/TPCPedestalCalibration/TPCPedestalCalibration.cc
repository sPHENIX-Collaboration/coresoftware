#include "TPCPedestalCalibration.h"

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

#include <TMath.h>
#include <TFile.h>
#include <TTree.h>

#include <memory>
#include <cassert>
#include <iostream>

/*************************************************************/
/*                TPC Pedestal Calibration                   */
/*               Thomas Marshall,Aditya Dash                 */
/*        rosstom@ucla.edu,aditya55@physics.ucla.edu         */
/*************************************************************/

TPCPedestalCalibration::TPCPedestalCalibration(const std::string &name)
 :SubsysReco("TPCPedestalCalibration")
 , m_fname(name)
{
  // reserve memory for max ADC samples
  m_adcSamples.resize(1024, 0);

  for(int fee_no=0;fee_no<26;fee_no++)
  {
    for(int channel_no=0;channel_no<256;channel_no++)
    {
      m_aveADCFeeChannel[fee_no][channel_no]=0.0;
      m_stdADCFeeChannel[fee_no][channel_no]=0.0;
      m_countsADCFeeChannel[fee_no][channel_no]=0.0;
      m_aliveArrayFeeChannel[fee_no][channel_no]=1;
    }
  }

}

int TPCPedestalCalibration::InitRun(PHCompositeNode *)
{
  m_file = TFile::Open(m_fname.c_str(), "recreate");
  assert(m_file->IsOpen());
 
  m_pedestalTree = new TTree("pedestalTree", "Each entry is one TPC Channel");

  m_pedestalTree->Branch("isAlive",&m_isAlive,"isAlive/I");
  m_pedestalTree->Branch("pedMean",&m_pedMean,"pedMean/F");
  m_pedestalTree->Branch("pedStd",&m_pedStd,"pedStd/F");
  m_pedestalTree->Branch("sector",&m_sector,"sector/I");
  m_pedestalTree->Branch("fee",&m_outFEE,"fee/I");
  m_pedestalTree->Branch("channel",&m_chan,"channel/I");
  m_pedestalTree->Branch("module",&m_module,"module/I");
  m_pedestalTree->Branch("slot",&m_slot,"slot/I");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TPCPedestalCalibration::process_event(PHCompositeNode *topNode)
{
  Event *_event = findNode::getClass<Event>(topNode, "PRDF");
  if (_event == nullptr)
  {
    std::cout << "TPCRawDataTree::Process_Event - Event not found" << std::endl;
    return -1;
  }
  if (_event->getEvtType() >= 8)  /// special events
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  for (int packet : m_packets)
  {
    if (Verbosity())
    {
      std::cout << __PRETTY_FUNCTION__ << " : decoding packet " << packet << std::endl;
    }

    m_packet = packet;

    std::unique_ptr<Packet> p (_event->getPacket(m_packet));
    if (!p)
    {
      if (Verbosity())
      {
        std::cout << __PRETTY_FUNCTION__ << " : missing packet " << packet << std::endl;
      }

      continue;
    }

    m_nWaveormInFrame = p->iValue(0, "NR_WF");
  
    for (int wf = 0; wf < m_nWaveormInFrame; wf++)
    {
      m_nSamples = p->iValue(wf, "SAMPLES");
      m_fee = p->iValue(wf, "FEE");
      m_Channel = p->iValue(wf, "CHANNEL");

      if(m_nSamples==0)
      {
        m_aliveArrayFeeChannel[m_fee][m_Channel]=0;
        continue;
      }

      assert(m_nSamples < (int) m_adcSamples.size());  // no need for movements in memory allocation
      for (int s = 0; s < m_nSamples; s++)
      {
        m_adcSamples[s] = p->iValue(wf, s);
      }

      bool dead = false;
      for(int adc_sam_no=0;adc_sam_no<m_nSamples;adc_sam_no++)
      {
        if (m_adcSamples[adc_sam_no] == 0 || TMath::IsNaN(float(m_adcSamples[adc_sam_no])))
        {
          m_aliveArrayFeeChannel[m_fee][m_Channel]=0;
        }
      }
      if (dead) {continue;}

      for(int adc_sam_no=0;adc_sam_no<m_nSamples;adc_sam_no++){
        m_aveADCFeeChannel[m_fee][m_Channel]+=m_adcSamples[adc_sam_no];
        m_stdADCFeeChannel[m_fee][m_Channel]+=pow(m_adcSamples[adc_sam_no],2);
        m_countsADCFeeChannel[m_fee][m_Channel]+=1.0;
      }
    }
  }  //   for (int packet : m_packets)

  return Fun4AllReturnCodes::EVENT_OK;
}

int TPCPedestalCalibration::EndRun(const int runnumber)
{
  std::cout << "TPCPedestalCalibration::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
 
  for(int fee_no=0;fee_no<26;fee_no++)
  {
    for(int channel_no=0;channel_no<256;channel_no++)
    {
      if(m_countsADCFeeChannel[fee_no][channel_no] != 0.0)
      {
        float temp1 = m_aveADCFeeChannel[fee_no][channel_no]/m_countsADCFeeChannel[fee_no][channel_no];
        float temp2 = m_stdADCFeeChannel[fee_no][channel_no]/m_countsADCFeeChannel[fee_no][channel_no];
        m_aveADCFeeChannel[fee_no][channel_no] = temp1;
        m_stdADCFeeChannel[fee_no][channel_no] = temp2;
      }
      else
      {
        m_aveADCFeeChannel[fee_no][channel_no] = 0.0;
        m_stdADCFeeChannel[fee_no][channel_no] = 0.0;
        m_aliveArrayFeeChannel[fee_no][channel_no]=0;
      }

      if(m_aveADCFeeChannel[fee_no][channel_no] > 200 || m_aveADCFeeChannel[fee_no][channel_no] < 10)
      {
        m_aliveArrayFeeChannel[fee_no][channel_no]=0;
      }

      m_pedMean=m_aveADCFeeChannel[fee_no][channel_no];
      m_pedStd=sqrt(m_stdADCFeeChannel[fee_no][channel_no] - pow(m_aveADCFeeChannel[fee_no][channel_no],2));
      m_isAlive=m_aliveArrayFeeChannel[fee_no][channel_no];
      m_chan=channel_no;
      m_outFEE=fee_no;
      m_module=mod_arr[fee_no];
      m_slot=slot_arr[fee_no];
      m_pedestalTree->Fill();
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TPCPedestalCalibration::End(PHCompositeNode *topNode)
{
  std::cout << "TPCPedestalCalibration::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  m_file->Write();

  std::cout << __PRETTY_FUNCTION__ << " : completed saving to " << m_file->GetName() << std::endl;
  m_file->ls();

  m_file->Close();  

  return Fun4AllReturnCodes::EVENT_OK;
}
