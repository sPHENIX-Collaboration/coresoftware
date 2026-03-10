#include "BcoLumiReco.h"

#include "BcoInfov1.h"

#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE
#include <phool/recoConsts.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>

#include <iostream>

BcoLumiReco::BcoLumiReco(const std::string &name)
  : SubsysReco(name)
{
  return;
}

int BcoLumiReco::Init(PHCompositeNode *topNode)
{
  int iret = CreateNodeTree(topNode);
  return iret;
}

int BcoLumiReco::InitRun(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
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
  BcoInfo *bcoinfo = findNode::getClass<BcoInfo>(topNode,"BCOINFO");
  if (!bcoinfo)
  {
    bcoinfo = new BcoInfov1();
    PHIODataNode<PHObject> *newnode = new PHIODataNode<PHObject>(bcoinfo,"BCOINFO","PHObject");
    dstNode->addNode(newnode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int BcoLumiReco::process_event(PHCompositeNode *topNode)
{
  static bool ifirst = true;
  Event* evt = findNode::getClass<Event>(topNode,"PRDF");
  if (evt)
  {
    evt->identify();
    if (evt->getEvtType() != DATAEVENT)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    Packet *packet = evt->getPacket(14001);
    uint64_t gtm_bco = packet->lValue(0, "BCO");
    std::cout << std::hex << "packet ival: 0x" << packet->lValue(0, "BCO")
	      << " uint64_t: 0x" << gtm_bco << std::dec << std::endl;
    push(gtm_bco);
    delete packet;
  }
  if (ifirst) // abort first event
  {
    ifirst = false;
    return Fun4AllReturnCodes::ABORTEVENT;
  }    
//  Fun4AllServer *se = Fun4AllServer::instance();
  SyncObject *syncobject = findNode::getClass<SyncObject>(topNode, syncdefs::SYNCNODENAME);
  if (!synccopy)
  {
    synccopy = dynamic_cast<SyncObject*> (syncobject->CloneMe()); // clone for second event
    tmpsync = dynamic_cast<SyncObject*> (synccopy->CloneMe()); // just to create this object
    return Fun4AllReturnCodes::ABORTEVENT; // and abort
  }
  BcoInfo *bcoinfo = findNode::getClass<BcoInfo>(topNode,"BCOINFO");

  std::cout << "current event is: " << syncobject->EventNumber() << "\n";
  std::cout << "saving as event: " << synccopy->EventNumber() << "\n";
  *tmpsync = *syncobject; // save current version
  *syncobject = *synccopy;
  *synccopy = *tmpsync;
  if (Verbosity() > 100)
  {
    std::cout << "current sync object\n";    
    syncobject->identify();
    std::cout << "next sync object\n";    
    synccopy->identify();
  }
  std::cout << std::hex;
  std::cout << "previous bco: " << get_previous_bco() << "\n";
  std::cout << "current  bco: " << get_current_bco() << "\n";
  std::cout << "future   bco: " << get_future_bco() << std::endl;
  std::cout << std::dec;
  bcoinfo->set_previous_bco(get_previous_bco());
  bcoinfo->set_current_bco(get_current_bco());
  bcoinfo->set_future_bco(get_future_bco());
  return Fun4AllReturnCodes::EVENT_OK;
}

void BcoLumiReco::push(uint64_t value)
{
  bco[0] = bco[1];
  bco[1] = bco[2];
  bco[2] = value;
}

  
