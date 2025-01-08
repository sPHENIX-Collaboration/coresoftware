#include "HcalCheck.h"

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
HcalCheck::HcalCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int HcalCheck::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HcalCheck::process_event(PHCompositeNode *topNode)
{
  CaloPacketContainer *hcalcont = findNode::getClass<CaloPacketContainer>(topNode, "HCALPackets");
  if (!hcalcont)
  {
    std::cout << "could not find HcalPacket node" << std::endl;
  }
  else
  {
    for (unsigned int i = 0; i < hcalcont->get_npackets(); i++)
    {
      if (ddump_enabled())
      {
	ddumppacket(hcalcont->getPacket(i));
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
