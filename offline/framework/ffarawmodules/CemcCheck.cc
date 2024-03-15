#include "CemcCheck.h"

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
CemcCheck::CemcCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int CemcCheck::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CemcCheck::process_event(PHCompositeNode *topNode)
{
  CaloPacketContainer *cemccont = findNode::getClass<CaloPacketContainer>(topNode, "CEMCPackets");
  if (!cemccont)
  {
    std::cout << "could not find CemcPacket node" << std::endl;
  }
  else
  {
    for (unsigned int i = 0; i < cemccont->get_npackets(); i++)
    {
      if (ddump_enabled())
      {
	ddumppacket(cemccont->getPacket(i));
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
