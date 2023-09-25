#include "InttCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/InttRawHitContainer.h>
#include <ffarawobjects/InttRawHit.h>

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
InttCheck::InttCheck(const std::string &name)
: SubsysReco(name)
{
}

//____________________________________________________________________________..
int InttCheck::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int InttCheck::process_event(PHCompositeNode *topNode)
{
  InttRawHitContainer *inttcont = findNode::getClass<InttRawHitContainer>(topNode,m_EvtNodeName);
  if (!inttcont)
  {
    std::cout << "could not find node " << m_EvtNodeName << std::endl;
  }
  else
  {
//  inttcont->identify();
    uint64_t refbco = std::numeric_limits<uint64_t>::max();
    bool ifirst = true;
    for (unsigned int i=0; i<inttcont->get_nhits(); i++)
    {
      InttRawHit *inh = inttcont->get_hit(i);
      if (ifirst)
      {
	refbco = inh->get_bco();
        ifirst = false;
      }
      else
      {
	if (refbco != inh->get_bco())
	{
	  std::cout << "scream, refbco: 0x" << std::hex << refbco
		    << " current bco: 0x" << inh->get_bco()
		    << std::dec << std::endl;
	}
      }
    }

  }
  return Fun4AllReturnCodes::EVENT_OK;
}
