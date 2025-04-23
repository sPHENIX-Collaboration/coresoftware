#include "SEpdCheck.h"

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
SEpdCheck::SEpdCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int SEpdCheck::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SEpdCheck::process_event(PHCompositeNode *topNode)
{
  CaloPacketContainer *sepdcont = findNode::getClass<CaloPacketContainer>(topNode, "SEPDPackets");
  if (!sepdcont)
  {
    std::cout << "could not find SEpdPacket node" << std::endl;
  }
  else
  {
    for (unsigned int i = 0; i < sepdcont->get_npackets(); i++)
    {
      if (ddump_enabled())
      {
        ddumppacket(sepdcont->getPacket(i));
      }
    }
    std::cout << "SEPD Evt no: " << sepdcont->getEvtSequence() << std::endl;
    for (unsigned int i = 0; i < sepdcont->get_npackets(); i++)
    {
      std::cout << "Packet " << sepdcont->getPacket(i)->getIdentifier()
                << " Evt no: " << sepdcont->getPacket(i)->getEvtSequence() << std::endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
