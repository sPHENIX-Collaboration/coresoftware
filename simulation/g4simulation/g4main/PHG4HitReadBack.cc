#include "PHG4HitReadBack.h"

#include "PHG4Hit.h"
#include "PHG4HitContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

class PHCompositeNode;

PHG4HitReadBack::PHG4HitReadBack(const std::string &name)
  : SubsysReco(name)
{
  return;
}

int PHG4HitReadBack::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *phc = findNode::getClass<PHG4HitContainer>(topNode, "PHG4Hit");
  if (!phc)
  {
    std::cout << "No PHG4Hit found" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  phc->identify();
  std::pair<PHG4HitContainer::ConstIterator, PHG4HitContainer::ConstIterator> hititer = phc->getHits(1);
  PHG4HitContainer::ConstIterator begin;
  PHG4HitContainer::ConstIterator end;
  PHG4HitContainer::ConstIterator it;
  begin = hititer.first;
  end = hititer.second;
  for (it = begin; it != end; ++it)
  {
    std::cout << "key: 0x" << std::hex << it->first << std::dec << std::endl;
    std::cout << "x: " << it->second->get_x(0) << std::endl;
  }
  std::cout << "detid: 2" << std::endl;
  hititer = phc->getHits(2);
  begin = hititer.first;
  end = hititer.second;
  for (it = begin; it != end; ++it)
  {
    std::cout << "key: 0x" << std::hex << it->first << std::dec << std::endl;
    std::cout << "x: " << it->second->get_x(0) << std::endl;
  }
  std::cout << "detid: 3" << std::endl;
  hititer = phc->getHits(3);
  begin = hititer.first;
  end = hititer.second;
  for (it = begin; it != end; ++it)
  {
    std::cout << "key: 0x" << std::hex << it->first << std::dec << std::endl;
    std::cout << "x: " << it->second->get_x(0) << std::endl;
  }
  //   phc = findNode::getClass<PHG4HitContainer>(topNode,"PHG4Hit2");
  //   if (!phc)
  //     {
  //       std::cout << "No PHG4Hit found" << std::endl;
  //       return ABORTEVENT;
  //     }
  //   phc->identify();
  return Fun4AllReturnCodes::EVENT_OK;
}
