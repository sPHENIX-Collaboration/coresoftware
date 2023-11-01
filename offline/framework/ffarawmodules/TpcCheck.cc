#include "TpcCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/TpcRawHitContainer.h>
#include <ffarawobjects/TpcRawHit.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <Event/packet.h>

#include <TSystem.h>

#include <iostream>  // for operator<<, endl, basic_ost...
#include <set>
#include <utility>   // for pair
#include <vector>    // for vector

//____________________________________________________________________________..
TpcCheck::TpcCheck(const std::string &name)
: SubsysReco(name)
{
}

//____________________________________________________________________________..
int TpcCheck::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TpcCheck::process_event(PHCompositeNode *topNode)
{
  TpcRawHitContainer *tpccont = findNode::getClass<TpcRawHitContainer>(topNode,m_EvtNodeName);
  if (!tpccont)
  {
    std::cout << "could not find node " << m_EvtNodeName << std::endl;
  }
  else
  {
  tpccont->identify();
  std::set<uint64_t> refbco;
    bool ifirst = true;
    for (unsigned int i=0; i<tpccont->get_nhits(); i++)
    {
      TpcRawHit *inh = tpccont->get_hit(i);
      if (ifirst)
      {
	refbco.insert(inh->get_bco());
	if (Verbosity() > 0)
	{
	  std::cout << "current bco: 0x" << std::hex << *refbco.begin()
		    << std::dec << std::endl;
	}
	if (bclk_seen.find(*refbco.begin()) != bclk_seen.end())
	{
	  std::cout << "bco 0x" << std::hex << *refbco.begin() << std::dec
		    << " seen before" << std::endl;
	}
	bclk_seen.insert(*refbco.begin());
        ifirst = false;
      }
      else
      {
	if (refbco.find(inh->get_bco()) == refbco.end())
	{
	  if (*refbco.begin() <= inh->get_bco() + bcorange)
	  {
	    refbco.insert(inh->get_bco());
	  }
	  else
	  {
	    std::cout << "scream, refbco: 0x" << std::hex << *refbco.begin()
		    << " current bco: 0x" << inh->get_bco()
		    << std::dec
		    << ", diff: " << inh->get_bco() - *refbco.begin()
		    << ", allowed diff: " << bcorange
<< std::endl;
	  }
	}
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
