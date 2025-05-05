#include "LL1Check.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/LL1Packet.h>
#include <ffarawobjects/LL1PacketContainer.h>

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
LL1Check::LL1Check(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int LL1Check::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int LL1Check::process_event(PHCompositeNode *topNode)
{
  LL1PacketContainer *ll1cont = findNode::getClass<LL1PacketContainer>(topNode, "LL1Packets");
  if (!ll1cont)
  {
    std::cout << "could not find LL1Packet node" << std::endl;
  }
  else
  {
    for (unsigned int i = 0; i < ll1cont->get_npackets(); i++)
    {
      if (ddump_enabled())
      {
        ddumppacket(ll1cont->getPacket(i));
      }
      //      ll1cont->getPacket(i)->identify();
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
