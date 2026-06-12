#include "StreamingBcoLumiReco.h"

#include "BcoInfo.h"
#include "StreamingBcoInfo.h"
#include "StreamingBcoInfov1.h"
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

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>  // for Packet

#include <TH1.h>

#include <iostream>

StreamingBcoLumiReco::StreamingBcoLumiReco(const std::string &name)
  : SubsysReco(name)
{
  hm = new Fun4AllHistoManager("bco_histos");
  Fun4AllServer *se = Fun4AllServer::instance();
  se->registerHistoManager(hm); 
  return;
}

int StreamingBcoLumiReco::Init(PHCompositeNode *topNode)
{
  int iret = CreateNodeTree(topNode);
  h_bco_diff = new TH1I("h_bco_diff", ";bco diff;", 3500, 0, 3500);
  std::string hist_name = "h_bco_diff_bit";
  for (int bit=0; bit<trigbits; ++bit)
  {
    hist_name += std::to_string(bit);
    h_bco_diff_trigbits[bit] = new TH1I(hist_name.c_str(), ";bco diff;", 3500, 0, 3500);
    hm->registerHisto(h_bco_diff_trigbits[bit]);
  }
  h_bco_tag = new TH1I("h_bco_tag", ";usable bco tag;", 2, -0.5, 1.5);
  hm->registerHisto(h_bco_diff);
  hm->registerHisto(h_bco_tag);

  return iret;
}

int StreamingBcoLumiReco::InitRun(PHCompositeNode * topNode)
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

int StreamingBcoLumiReco::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " DST Node is missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  StreamingBcoInfo *streaming_bco_info = findNode::getClass<StreamingBcoInfo>(topNode, "STREAMINGBCOINFO");
  if (!streaming_bco_info)
  {
    streaming_bco_info = new StreamingBcoInfov1();
    PHIODataNode<PHObject> *bconode = new PHIODataNode<PHObject>(streaming_bco_info, "STREAMINGBCOINFO", "PHObject");
    dstNode->addNode(bconode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int StreamingBcoLumiReco::process_event(PHCompositeNode *topNode)
{
  BcoInfo *bcoinfo = findNode::getClass<BcoInfo>(topNode, "BCOINFO");
  SyncObject *syncobject = findNode::getClass<SyncObject>(topNode, syncdefs::SYNCNODENAME);
  //Gl1Packet *gl1packet = findNode::getClass<Gl1Packet>(topNode, 14001);
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
    uint64_t gl1_scaledvec = packet->lValue(0, "ScaledVector");
    //uint64_t gl1_livevec = packet->lValue(0, "TriggerVector");   

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

    if (Verbosity() > 2)
    {
      if (!syncobject)
      {
        std::cout << PHWHERE << " SyncObject missing" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      std::cout << "Event No: " << syncobject->EventNumber() /*<< std::hex*/
                << " gl1  bco: " << gtm_bco <<std::dec << std::endl;
    }
    if (bcoinfo)
    {
      if (Verbosity() > 2)
      {
        std::cout << "prev event: " << bcoinfo->get_previous_evtno() /*<< std::hex*/
                  << " bco: " << bcoinfo->get_previous_bco() << std::dec << std::endl;
        std::cout << "curr event: " << bcoinfo->get_current_evtno() /*<< std::hex*/
                  << " bco: " << bcoinfo->get_current_bco() << std::dec << std::endl;
        std::cout << "futu event: " << bcoinfo->get_future_evtno() /*<< std::hex*/
                  << " bco: " << bcoinfo->get_future_bco() << std::dec << std::endl;
      }

      StreamingBcoInfo *streaming_bco_info = findNode::getClass<StreamingBcoInfo>(topNode, "STREAMINGBCOINFO");
      if (!streaming_bco_info)
      {
        std::cout << PHWHERE << " STREAMINGBCOINFO node missing" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      m_bco = bcoinfo->get_current_bco();
      if (gtm_bco != m_bco) { std::cout << "BCO MISMATCH!!! : gtm_bco : " << gtm_bco << " m_bco " << m_bco << std::endl;}
      uint64_t bco_prev = bcoinfo->get_previous_bco();
      uint64_t bco_futu = bcoinfo->get_future_bco();
      uint64_t bco_diff_prev = m_bco - bco_prev;
      uint64_t bco_diff_futu = bco_futu - m_bco;

      // special case if BCO is within 20 of previous BCO?
      if (bco_diff_prev < m_default_positive_window_length)
      {
        m_usable_bco_tag = true;
      }
      else
      {
        m_usable_bco_tag = false;
      }
      if (bco_diff_futu < m_default_positive_window_length)
      {
        // double check boundaries for overlap!!
        m_bco_streaming_window = std::make_pair(get_bco() - m_default_negative_window_length, bco_futu - m_default_negative_window_length + 1);
      }
      else
      {
        m_bco_streaming_window = std::make_pair(get_bco() - m_default_negative_window_length, get_bco() + m_default_positive_window_length);
      }
      if (Verbosity() > 2)
      {
        std::cout << "bco_diff_prev : " << bco_diff_prev << std::endl;
        std::cout << "bco_diff_futu : " << bco_diff_futu << std::endl; 
      }
      h_bco_diff->Fill(bco_diff_prev);
      h_bco_tag->Fill(m_usable_bco_tag);
      for (int bit=0; bit<trigbits; ++bit)
      {
        bool trigger_fired = ((gl1_scaledvec >> static_cast<uint64_t>(bit)) & 0x1U) == 0x1U;
        //bool scaled_trigger_fired = ((gl1_scaledvec >> bit) & 0x1U) == 0x1U;

        if (trigger_fired) 
        {
          h_bco_diff_trigbits[bit]->Fill(bco_diff_prev);
        }
      }
      // Double check Zhiwan's logic for assigning the adjusted bunch!
      int lower = m_bco_streaming_window.first - m_bco;
      int upper = m_bco_streaming_window.second - m_bco;
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
          if(i!=0 || m_usable_bco_tag) 
          {
            m_bunchnumber_crossings[adjusted_bunch] += 1;
          }
          //else if (m_usable_bco_tag)
          //{
          //  m_bunchnumber_crossings[adjusted_bunch] += 1;
          //}
      }

      streaming_bco_info->set_bco(get_bco());
      streaming_bco_info->set_usable_bco_tag(get_usable_bco_tag());
      streaming_bco_info->set_bco_streaming_window(get_bco_streaming_window());
      if (syncobject)
      {
        streaming_bco_info->set_evtno(syncobject->EventNumber());
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int StreamingBcoLumiReco::EndRun(int /*runnumber*/)
{
  uint64_t rawgl1scalers_per_bunch = m_rawgl1scaler/120.;
  for (int i=0; i<m_bunches; i++)
  {
    m_bunchnumber_lumi_raw[i] = ( m_bunchnumber_crossings[i] * m_bunchnumber_MBDNS_raw[i] ) / ( rawgl1scalers_per_bunch * m_xsec_MBDNS );
    m_bunchnumber_lumi_live[i] = ( m_bunchnumber_crossings[i] * m_bunchnumber_MBDNS_live[i]) / ( rawgl1scalers_per_bunch * m_xsec_MBDNS );
    m_bunchnumber_lumi_scaled[i] = ( m_bunchnumber_crossings[i] * m_bunchnumber_MBDNS_scaled[i] ) / ( rawgl1scalers_per_bunch * m_xsec_MBDNS );

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
  m_streaming_lumi_info->set_lumi_raw(get_lumi_raw());
  m_streaming_lumi_info->set_lumi_live(get_lumi_live());
  m_streaming_lumi_info->set_lumi_scaled(get_lumi_scaled());
  if (Verbosity() > 1)
  {
    std::cout << "MBD xsec : " << m_xsec_MBDNS << std::endl;
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
