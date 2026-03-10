#include "BcoLumiCheck.h"

#include "BcoInfo.h"

#include <ffaobjects/SyncDefs.h>
#include <ffaobjects/SyncObject.h>

#include <ffarawobjects/Gl1Packet.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <iostream>

BcoLumiCheck::BcoLumiCheck(const std::string &name)
  : SubsysReco(name)
{
  return;
}

int BcoLumiCheck::Init(PHCompositeNode *topNode)
{
  int iret = CreateNodeTree(topNode);
  return iret;
}

int BcoLumiCheck::InitRun(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int BcoLumiCheck::CreateNodeTree(PHCompositeNode *topNode)
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

int BcoLumiCheck::process_event(PHCompositeNode *topNode)
{
  BcoInfo *bcoinfo = findNode::getClass<BcoInfo>(topNode, "BCOINFO");
  SyncObject *syncobject = findNode::getClass<SyncObject>(topNode, syncdefs::SYNCNODENAME);
  Gl1Packet *gl1packet = findNode::getClass<Gl1Packet>(topNode, 14001);
  if (gl1packet)
  {
    std::cout << "Event No: " << syncobject->EventNumber() << std::hex
              << " gl1:  bco 0x" << gl1packet->lValue(0, "BCO") << std::dec << std::endl;
    if (bcoinfo)
    {
      std::cout << "prev event: " << bcoinfo->get_previous_evtno() << std::hex
                << " bco: 0x" << bcoinfo->get_previous_bco() << std::dec << std::endl;
      std::cout << "curr event: " << bcoinfo->get_current_evtno() << std::hex
                << " bco: 0x" << bcoinfo->get_current_bco() << std::dec << std::endl;
      std::cout << "futu event: " << bcoinfo->get_future_evtno() << std::hex
                << " bco: 0x" << bcoinfo->get_future_bco() << std::dec << std::endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
