#include "EvtCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <fun4allraw/PacketMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <Event/packet.h>

#include <TSystem.h>

#include <iostream>  // for operator<<, endl, basic_ost...
#include <utility>   // for pair
#include <vector>    // for vector

//____________________________________________________________________________..
EvtCheck::EvtCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int EvtCheck::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EvtCheck::process_event(PHCompositeNode *topNode)
{
  PacketMap *pktmap = findNode::getClass<PacketMap>(topNode,m_EvtNodeName);
  if (!pktmap)
  {
    std::cout << "could not find node " << m_EvtNodeName << std::endl;
  }
  else
  {
  pktmap->identify();
  }
  PacketMap::PacketRange pktrange = pktmap->first_last_packet();
  for (auto iter = pktrange.first; iter != pktrange.second; iter++)
  {
    std::cout << "Packet " << iter->first << std::endl;
    for (auto bclkiter : iter->second.m_BeamClockSet)
    {
      std::cout << "bclks: 0x" << std::hex << bclkiter << std::dec << std::endl;
    }
    for (auto pktiter : iter->second.m_PacketVector)
    {
      std::cout << "Packet " <<  pktiter->getIdentifier() << " at "
		<< std::hex << pktiter << std::dec << std::endl;
    }

  }
  return Fun4AllReturnCodes::EVENT_OK;
}

