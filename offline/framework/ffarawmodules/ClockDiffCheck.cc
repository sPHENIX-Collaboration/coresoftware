#include "ClockDiffCheck.h"

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

#include <TH1.h>
#include <TSystem.h>

#include <iostream>  // for operator<<, endl, basic_ost...
#include <utility>   // for pair
#include <vector>    // for vector

//____________________________________________________________________________..
ClockDiffCheck::ClockDiffCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int ClockDiffCheck::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ClockDiffCheck::process_event(PHCompositeNode *topNode)
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
      unsigned int packetid = cemccont->getPacket(i)->getIdentifier();
      unsigned int clk = cemccont->getPacket(i)->getBCO();
	if (m_PacketStuffMap.find(packetid) == m_PacketStuffMap.end())
	{
	  std::string hname = "clkdiff" + std::to_string(packetid);
	  TH1 *h1 = new TH1F(hname.c_str(), hname.c_str(),100,0,99);
	  m_PacketStuffMap[packetid] = std::make_tuple(clk, h1);
	}
      if (ddump_enabled())
      {
	ddumppacket(cemccont->getPacket(i));
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
