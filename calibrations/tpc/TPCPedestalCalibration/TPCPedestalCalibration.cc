#include "TPCPedestalCalibration.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/recoConsts.h>
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
 , m_writeToCDB(false)
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
  m_cdbttree = new CDBTTree(m_fname);  

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
      if (m_firstBCO == true)
      {
        m_BCO = p->iValue(wf, "BCO");
        m_firstBCO = false;
      }

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
      
      m_cdbttree->SetIntValue(fee_no*256 + channel_no,"isAlive",m_isAlive);
      m_cdbttree->SetFloatValue(fee_no*256 + channel_no,"pedMean",m_pedMean);
      m_cdbttree->SetFloatValue(fee_no*256 + channel_no,"pedStd",m_pedStd);
      m_cdbttree->SetIntValue(fee_no*256 + channel_no,"sector",m_sector);
      m_cdbttree->SetIntValue(fee_no*256 + channel_no,"fee",m_outFEE);
      m_cdbttree->SetIntValue(fee_no*256 + channel_no,"channel",m_chan);
      m_cdbttree->SetIntValue(fee_no*256 + channel_no,"module",m_module);
      m_cdbttree->SetIntValue(fee_no*256 + channel_no,"slot",m_slot);
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void TPCPedestalCalibration::CDBInsert()
{
   recoConsts *rc = recoConsts::instance();
   
   CDBUtils *cdbInsert = new CDBUtils(rc->get_StringFlag("CDB_GLOBALTAG"));
   cdbInsert->createPayloadType("TPCPedestalCalibration");
   cdbInsert->insertPayload("TPCPedestalCalibration",m_fname,m_BCO); // uses first BCO value from first waveform
   
   return;
}

//____________________________________________________________________________..
int TPCPedestalCalibration::End(PHCompositeNode *topNode)
{
  m_cdbttree->Commit();
  if (Verbosity())
  {
    m_cdbttree->Print();
  }
  m_cdbttree->WriteCDBTTree();
  delete m_cdbttree;
 
  if (m_writeToCDB)
  {
    std::cout << "Inserting file " << m_fname << " in CDB under username " << m_username << std::endl;
    CDBInsert(); 
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
