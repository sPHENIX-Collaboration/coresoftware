#include "Gl1Check.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/Gl1Packet.h>

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
Gl1Check::Gl1Check(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int Gl1Check::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int Gl1Check::process_event(PHCompositeNode *topNode)
{
  Gl1Packet *gl1cont = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  if (!gl1cont)
  {
    std::cout << "could not find Gl1Packet node" << std::endl;
  }
  else
  {
    gl1cont->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
