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
  inttcont->identify();

  for (unsigned int i=0; i<inttcont->get_nhits(); i++)
  {
    inttcont->get_hit(i)->identify();
  }

  }
  return Fun4AllReturnCodes::EVENT_OK;
}

