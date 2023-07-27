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
  PacketMap::PacketListRange pktrange = pktmap->first_last_packet();
  for (auto iter = pktrange.first; iter != pktrange.second; iter++)
  {
    std::cout << "Packet " << iter->first << std::endl;
    for (auto bclkiter : iter->second.m_BeamClockSet)
    {
      std::cout << "bclks: 0x" << std::hex << bclkiter << std::dec << std::endl;
      for (auto pktiter : iter->second.m_PacketVector)
      {
	std::cout << "Packet " <<  pktiter->getIdentifier() << " at "
		  << std::hex << pktiter << std::dec << std::endl;
	int num_hits = pktiter->iValue(0, "NR_HITS");
        for (int j = 0; j < num_hits; j++)
        {
          int FEE = pktiter->iValue(j, "FEE");
          uint64_t gtm_bco = pktiter->lValue(j, "BCO");
	  if ((bclkiter & 0xFFFFFFFFFF) == gtm_bco)
	  {
	    std::cout << "FEE " << FEE << " bclk 0x" << std::hex
		      << gtm_bco << std::dec << std::endl;
	  }
	}
      }
    }

  }
  return Fun4AllReturnCodes::EVENT_OK;
}

