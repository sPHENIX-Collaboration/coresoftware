#include "EventNumberCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <Event/Event.h>
#include <Event/oncsEvent.h>

#include <TSystem.h>

#include <iostream>  // for operator<<, endl, basic_ost...
#include <utility>   // for pair
#include <vector>    // for vector

//____________________________________________________________________________..
EventNumberCheck::EventNumberCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int EventNumberCheck::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int EventNumberCheck::process_event(PHCompositeNode *topNode)
{
    Event *evt = findNode::getClass<Event>(topNode,m_MyPrdfNode);
    evt->identify();
  int nw = evt->getPacketList(plist, 10000);
  if (nw >= 10000)
  {
    std::cout << "Packet array too small, need " << nw << " entries" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  for (int i=0; i<nw; i++)
  {
    std::cout << "packet " << plist[i]->getIdentifier() << ", evt nr "
	 <<   plist[i]->iValue(0, "EVTNR") << std::endl;
    delete plist[i];
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

