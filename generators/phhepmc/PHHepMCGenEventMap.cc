#include "PHHepMCGenEventMap.h"

#include "PHHepMCGenEvent.h"

#include <TSystem.h>

#include <cassert>
#include <cstdlib>           // for exit
#include <iterator>           // for reverse_iterator
#include <utility>            // for pair, make_pair

PHHepMCGenEventMap::PHHepMCGenEventMap(const PHHepMCGenEventMap& eventmap)
{
  for( const auto& pair:eventmap.get_map() )
  { _map.insert(std::make_pair(pair.first, static_cast<PHHepMCGenEvent*>(pair.second->CloneMe()))); }
}

PHHepMCGenEventMap& PHHepMCGenEventMap::operator=(const PHHepMCGenEventMap& eventmap)
{
  Reset();
  for( const auto& pair:eventmap.get_map() )
  { _map.insert(std::make_pair(pair.first, static_cast<PHHepMCGenEvent*>(pair.second->CloneMe()))); }
  return *this;
}

PHHepMCGenEventMap::~PHHepMCGenEventMap()
{
  Reset();
}

void PHHepMCGenEventMap::Reset()
{
  // delete all events
  for( const auto& pair:_map )
  {
    delete pair.second;
  }
  
  // clear map
  _map.clear();
}

void PHHepMCGenEventMap::identify(std::ostream& os) const
{
  os << "PHHepMCGenEventMap: size = " << _map.size() << std::endl;

  for (const auto& evt : _map)
  {
    std::cout << "Event[" << evt.first << "] : ";
    assert(evt.second);
    evt.second->identify();
  }

  return;
}

const PHHepMCGenEvent* PHHepMCGenEventMap::get(int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

PHHepMCGenEvent* PHHepMCGenEventMap::get(int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end()) return nullptr;
  return iter->second;
}

PHHepMCGenEvent* PHHepMCGenEventMap::insert_active_event(const PHHepMCGenEvent* event)
{
  const int index = _map.empty() ? 1 : (_map.rbegin()->first+1);
  auto newEvent = event ? static_cast<PHHepMCGenEvent*>(event->CloneMe()): new PHHepMCGenEvent();
  newEvent->set_embedding_id(index);
  _map.insert(std::make_pair(index,newEvent));
  return newEvent;
}

PHHepMCGenEvent* PHHepMCGenEventMap::insert_background_event(const PHHepMCGenEvent* event)
{
  const int index = _map.empty() ? -1 : (_map.begin()->first-1);
  auto newEvent = event ? static_cast<PHHepMCGenEvent*>(event->CloneMe()): new PHHepMCGenEvent();
  newEvent->set_embedding_id(index);
  _map.insert(std::make_pair(index,newEvent));
  return newEvent;
}

PHHepMCGenEvent* PHHepMCGenEventMap::insert_event(const int index, const PHHepMCGenEvent* event)
{
  auto newEvent = event ? static_cast<PHHepMCGenEvent*>(event->CloneMe()): new PHHepMCGenEvent();
  newEvent->set_embedding_id(index);
  auto ret =  _map.insert(std::make_pair(index, newEvent ));
  if (ret.second == false)
  {
    std::cout 
      << "PHHepMCGenEventMap::insert_event - Fatal Error -"
      << "embedding ID " << index << " is already used in the PHHepMCGenEventMap. Print map:";
    identify();
    gSystem->Exit(10);
  }
  return newEvent;
}
