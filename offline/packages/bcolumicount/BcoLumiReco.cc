#include "BcoLumiReco.h"

#include "BcoInfo.h"
#include "BcoInfov1.h"

#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>  // for Packet
#
#include <iostream>

BcoLumiReco::BcoLumiReco(const std::string &name)
  : SubsysReco(name)
{
  return;
}

BcoLumiReco::~BcoLumiReco()
{
  delete m_synccopy;
  delete m_tmpsync;
}

int BcoLumiReco::Init(PHCompositeNode *topNode)
{
  int iret = CreateNodeTree(topNode);
  return iret;
}

int BcoLumiReco::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " DST Node is missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  BcoInfo *bcoinfo = findNode::getClass<BcoInfo>(topNode, "BCOINFO");
  if (!bcoinfo)
  {
    bcoinfo = new BcoInfov1();
    PHIODataNode<PHObject> *newnode = new PHIODataNode<PHObject>(bcoinfo, "BCOINFO", "PHObject");
    dstNode->addNode(newnode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int BcoLumiReco::process_event(PHCompositeNode *topNode)
{
  static bool ifirst = true;
  SyncObject *syncobject = findNode::getClass<SyncObject>(topNode, syncdefs::SYNCNODENAME);
  Event *evt = findNode::getClass<Event>(topNode, "PRDF");
  if (evt)
  {
    if (Verbosity() > 1)
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
    if (Verbosity() > 1)
    {
      std::cout << std::hex << "packet ival: 0x" << packet->lValue(0, "BCO")
                << " uint64_t: 0x" << gtm_bco << std::dec << std::endl;
    }
    push_bco(gtm_bco);
    delete packet;
  }
  if (syncobject)
  {
    push_evtno(syncobject->EventNumber());
  }
  if (ifirst)  // abort first event since it does not have a previous bco
  {
    ifirst = false;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!m_synccopy)
  {
    m_synccopy = dynamic_cast<SyncObject *>(syncobject->CloneMe());  // clone for second event
    m_tmpsync = dynamic_cast<SyncObject *>(m_synccopy->CloneMe());   // just to create this object
    return Fun4AllReturnCodes::ABORTEVENT;                           // and abort
  }
  BcoInfo *bcoinfo = findNode::getClass<BcoInfo>(topNode, "BCOINFO");

  if (Verbosity() > 0)
  {
    std::cout << "current event is: " << syncobject->EventNumber() << "\n";
    std::cout << "saving as event: " << m_synccopy->EventNumber() << "\n";
  }
  // here we store the current sync object and overwrite its content with the cached copy
  *m_tmpsync = *syncobject;   // save current version in tmp
  *syncobject = *m_synccopy;  // copy previously cached version
  *m_synccopy = *m_tmpsync;   // cache current version
  if (Verbosity() > 0)
  {
    std::cout << std::hex;
    std::cout << "previous bco: " << get_previous_bco() << "\n";
    std::cout << "current  bco: " << get_current_bco() << "\n";
    std::cout << "future   bco: " << get_future_bco() << std::endl;
    std::cout << std::dec;
  }
  bcoinfo->set_previous_bco(get_previous_bco());
  bcoinfo->set_current_bco(get_current_bco());
  bcoinfo->set_future_bco(get_future_bco());
  bcoinfo->set_previous_evtno(get_previous_evtno());
  bcoinfo->set_current_evtno(get_current_evtno());
  bcoinfo->set_future_evtno(get_future_evtno());
  return Fun4AllReturnCodes::EVENT_OK;
}

void BcoLumiReco::push_bco(uint64_t value)
{
  m_bco[0] = m_bco[1];
  m_bco[1] = m_bco[2];
  m_bco[2] = value;
}

void BcoLumiReco::push_evtno(int value)
{
  m_evtno[0] = m_evtno[1];
  m_evtno[1] = m_evtno[2];
  m_evtno[2] = value;
}
