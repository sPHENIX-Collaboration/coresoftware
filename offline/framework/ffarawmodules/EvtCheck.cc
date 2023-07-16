#include "EvtCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <fun4allraw/PacketList.h>

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
  PacketList *pktlist = findNode::getClass<PacketList>(topNode,m_EvtNodeName);
  if (!pktlist)
  {
    std::cout << "could not find node " << m_EvtNodeName << std::endl;
  }
  else
  {
  pktlist->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

