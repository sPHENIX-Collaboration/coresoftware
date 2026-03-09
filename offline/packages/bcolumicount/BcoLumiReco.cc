#include "BcoLumiReco.h"

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

  return Fun4AllReturnCodes::EVENT_OK;
}

int BcoLumiReco::process_event(PHCompositeNode *topNode)
{
  static bool ifirst = true;
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
  Event* evt = findNode::getClass<Event>(topNode,"PRDF");
  if (evt)
  {
    evt->identify();
      Packet *packet = evt->getPacket(14001);
uint64_t gtm_bco = packet->lValue(0, "BCO");
push(gtm_bco);
delete packet;
  }
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
  return Fun4AllReturnCodes::EVENT_OK;
}

void BcoLumiReco::push(uint64_t value)
{
  bco[0] = bco[1];
  bco[1] = bco[2];
  bco[2] = value;
}

  
