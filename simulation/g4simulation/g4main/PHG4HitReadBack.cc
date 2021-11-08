#include "PHG4HitReadBack.h"

#include "PHG4Hit.h"
#include "PHG4HitContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>

#include <iostream>
#include <map>                           // for _Rb_tree_const_iterator
#include <utility>                       // for pair

class PHCompositeNode;

using namespace std;

PHG4HitReadBack::PHG4HitReadBack(const string &name): SubsysReco(name)
{
  return;
}

int
PHG4HitReadBack::InitRun(PHCompositeNode */*topNode*/)
{
  return 0;
}

int
PHG4HitReadBack::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *phc = findNode::getClass<PHG4HitContainer>(topNode,"PHG4Hit");
  if (!phc)
    {
      cout << "No PHG4Hit found" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  phc->identify();
  std::pair<PHG4HitContainer::ConstIterator,PHG4HitContainer::ConstIterator> hititer = phc->getHits(1);
  PHG4HitContainer::ConstIterator begin,end,it;
  begin = hititer.first;
  end = hititer.second;
  for (it=begin; it != end; ++it)
    {
      cout << "key: 0x" << hex << it->first << dec << endl;
      cout << "x: " << it->second->get_x(0) << endl;
    }
  cout << "detid: 2" << endl;
  hititer = phc->getHits(2);
  begin = hititer.first;
  end = hititer.second;
  for (it=begin; it != end; ++it)
    {
      cout << "key: 0x" << hex << it->first << dec << endl;
      cout << "x: " << it->second->get_x(0) << endl;
    }
  cout << "detid: 3" << endl;
  hititer = phc->getHits(3);
  begin = hititer.first;
  end = hititer.second;
  for (it=begin; it != end; ++it)
    {
      cout << "key: 0x" << hex << it->first << dec << endl;
      cout << "x: " << it->second->get_x(0) << endl;
    }
  //   phc = findNode::getClass<PHG4HitContainer>(topNode,"PHG4Hit2");
  //   if (!phc)
  //     {
  //       cout << "No PHG4Hit found" << endl;
  //       return ABORTEVENT;
  //     }
  //   phc->identify();
  return Fun4AllReturnCodes::EVENT_OK;
}
