#include "MbdCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <TSystem.h>

#include <iostream>  // for operator<<, endl, basic_ost...
#include <utility>   // for pair
#include <vector>    // for vector

//____________________________________________________________________________..
MbdCheck::MbdCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int MbdCheck::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int MbdCheck::process_event(PHCompositeNode *topNode)
{
  CaloPacketContainer *mbdcont = findNode::getClass<CaloPacketContainer>(topNode, "MBDPackets");
  if (!mbdcont)
  {
    std::cout << "could not find MbdPacket node" << std::endl;
  }
  else
  {
    for (unsigned int i = 0; i < mbdcont->get_npackets(); i++)
    {
      if (ddump_enabled())
      {
	ddumppacket(mbdcont->getPacket(i));
      }
      mbdcont->getPacket(i)->identify();
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
