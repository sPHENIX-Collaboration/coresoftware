#include "StreamingLumiReco.h"

#include "StreamingBcoInfo.h"
#include "StreamingLumiInfo.h"
#include "StreamingLumiInfov1.h"


#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

#include <ffarawobjects/Gl1Packet.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllServer.h>


#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <phool/sphenix_constants.h> // for MDB_NS_xsec


#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>  // for Packet

#include <TH1.h>

#include <iostream>

StreamingLumiReco::StreamingLumiReco(const std::string &name)
  : SubsysReco(name)
{
  return;
}

int StreamingLumiReco::Init(PHCompositeNode *topNode)
{
  int iret = CreateNodeTree(topNode);
  return iret;
}

int StreamingLumiReco::InitRun(PHCompositeNode * topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode;
  runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << " Run Node is missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_streaming_lumi_info = findNode::getClass<StreamingLumiInfo>(topNode, "STREAMINGLUMIINFO");
  if (!m_streaming_lumi_info)
  {
    m_streaming_lumi_info = new StreamingLumiInfov1();
    PHIODataNode<PHObject> *luminode = new PHIODataNode<PHObject>(m_streaming_lumi_info, "STREAMINGLUMIINFO", "PHObject");
    runNode->addNode(luminode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int StreamingLumiReco::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " DST Node is missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int StreamingLumiReco::process_event(PHCompositeNode *topNode)
{
  StreamingBcoInfo *streaming_bcoinfo = findNode::getClass<StreamingBcoInfo>(topNode, "STREAMINGBCOINFO");
  Event *evt = findNode::getClass<Event>(topNode, "PRDF");
  if (evt)
  {
    if (Verbosity() > 2)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    Packet *packet = evt->getPacket(14001);
    if (!packet)
    {
      if (Verbosity() > 0)
      {
        std::cout << "no gl1 packet 14001" << std::endl;
        evt->identify();
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    uint64_t gtm_bco = packet->lValue(0, "BCO"); 

    int bunchno = packet->lValue(0,"BunchNumber");
    if (bunchno < 0 || bunchno >= m_bunches)
    {
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << " invalid bunch number: " << bunchno << std::endl;
      }
      delete packet;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    // SYNTAX TAKEN FROM ZHIWANS CODE, why = and not +=? If this is correct it seems like a waste to call it for every event (would just need it for the last event in a particular crossing?)
    m_bunchnumber_MBDNS_raw[bunchno] = packet->lValue(0, "GL1PRAW");
    m_bunchnumber_MBDNS_live[bunchno] = packet->lValue(0, "GL1PLIVE");
    m_bunchnumber_MBDNS_scaled[bunchno] = packet->lValue(0, "GL1PSCALED");


    if(packet->lValue(0, 0))
    {
      m_rawgl1scaler = packet->lValue(0, 0);
    }

    delete packet;

    if (streaming_bcoinfo)
    {
      if (gtm_bco != streaming_bcoinfo->get_bco()) { std::cout << "BCO MISMATCH!!! : gtm_bco : " << gtm_bco << " bco " << streaming_bcoinfo->get_bco() << std::endl;}

      // Double check Zhiwan's logic for assigning the adjusted bunch!
      int lower = streaming_bcoinfo->get_bco_streaming_window().first - streaming_bcoinfo->get_bco();
      int upper = streaming_bcoinfo->get_bco_streaming_window().second - streaming_bcoinfo->get_bco();
      for(int i = lower; i< upper;i++)
      {
        int adjusted_bunch = bunchno + i;
          while (adjusted_bunch < 0)
          {
              adjusted_bunch += 120;
          }
          while (adjusted_bunch > 119)
          {
              adjusted_bunch -= 120;
          }
          // ABORT GAP!
          if (adjusted_bunch>110) { continue; }

          // Make sure this is the correct way to count crossings! Need to zero out for each run!
          if(i!=0 || streaming_bcoinfo->get_usable_bco_tag())
          {
            m_bunchnumber_crossings[adjusted_bunch] += 1;
          }
          //else if (m_usable_bco_tag)
          //{
          //  m_bunchnumber_crossings[adjusted_bunch] += 1;
          //}
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int StreamingLumiReco::EndRun(int /*runnumber*/)
{
  uint64_t rawgl1scalers_per_bunch = m_rawgl1scaler/120.;
  for (int i=0; i<m_bunches; i++)
  {
    m_bunchnumber_lumi_raw[i] = ( m_bunchnumber_crossings[i] * m_bunchnumber_MBDNS_raw[i] ) / ( rawgl1scalers_per_bunch * sphenix_constants::m_xsec_MBDNS );
    m_bunchnumber_lumi_live[i] = ( m_bunchnumber_crossings[i] * m_bunchnumber_MBDNS_live[i]) / ( rawgl1scalers_per_bunch * sphenix_constants::m_xsec_MBDNS );
    m_bunchnumber_lumi_scaled[i] = ( m_bunchnumber_crossings[i] * m_bunchnumber_MBDNS_scaled[i] ) / ( rawgl1scalers_per_bunch * sphenix_constants::m_xsec_MBDNS );

    m_lumi_raw += m_bunchnumber_lumi_raw[i];
    m_lumi_live += m_bunchnumber_lumi_live[i];
    m_lumi_scaled += m_bunchnumber_lumi_scaled[i];

    if (Verbosity() > 1)
    {
     std::cout << "bunchno : " << i << " lumi_raw : " << m_bunchnumber_lumi_raw[i] << std::endl;
    }
  }
  if (!m_streaming_lumi_info)
  {
    std::cout << PHWHERE << " STREAMINGLUMIINFO node missing in EndRun" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_streaming_lumi_info->set_bunchnumber_lumi_raw(get_bunchnumber_lumi_raw());
  m_streaming_lumi_info->set_bunchnumber_lumi_live(get_bunchnumber_lumi_live());
  m_streaming_lumi_info->set_bunchnumber_lumi_scaled(get_bunchnumber_lumi_scaled());

  m_streaming_lumi_info->set_lumi_raw(get_lumi_raw());
  m_streaming_lumi_info->set_lumi_live(get_lumi_live());
  m_streaming_lumi_info->set_lumi_scaled(get_lumi_scaled());
  
  if (Verbosity() > 1)
  {
    std::cout << "MBD xsec : " << sphenix_constants::m_xsec_MBDNS << std::endl;
    std::cout << "total lumi (raw) : " << m_lumi_raw << std::endl;
  }

  m_bunchnumber_lumi_raw.fill(0);
  m_bunchnumber_lumi_live.fill(0);
  m_bunchnumber_lumi_scaled.fill(0);
  m_bunchnumber_crossings.fill(0);
  m_bunchnumber_MBDNS_raw.fill(0);
  m_bunchnumber_MBDNS_live.fill(0);
  m_bunchnumber_MBDNS_scaled.fill(0);
  m_lumi_raw = 0.;
  m_lumi_live = 0.;
  m_lumi_scaled = 0.;
  m_rawgl1scaler = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}
