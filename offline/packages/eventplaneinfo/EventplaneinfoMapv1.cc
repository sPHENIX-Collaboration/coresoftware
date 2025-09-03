#include "EventplaneinfoMapv1.h"

#include "Eventplaneinfo.h"
#include "EventplaneinfoMap.h"

#include <iterator>  // for reverse_iterator
#include <utility>   // for pair, make_pair

EventplaneinfoMapv1::~EventplaneinfoMapv1()
{
  EventplaneinfoMapv1::clear();
}

void EventplaneinfoMapv1::identify(std::ostream& os) const
{
  os << "EventplaneinfoMapv1: size = " << _map.size() << std::endl;
  for (const auto& m : _map)
  {
    m.second->identify(os);
  }
  return;
}

void EventplaneinfoMapv1::clear()
{
  for (const auto& iter : _map)
  {
    delete iter.second;
  }
  _map.clear();
  return;
}

const Eventplaneinfo* EventplaneinfoMapv1::get(unsigned int id) const
{
  ConstIter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

Eventplaneinfo* EventplaneinfoMapv1::get(unsigned int id)
{
  Iter iter = _map.find(id);
  if (iter == _map.end())
  {
    return nullptr;
  }
  return iter->second;
}

Eventplaneinfo* EventplaneinfoMapv1::insert(Eventplaneinfo* clus, const EventplaneinfoMap::EPTYPE id)
{
  auto [iter, inserted] = _map.insert(std::make_pair(id, clus));
  if (!inserted)
  {
    delete iter->second;
    iter->second = clus;
  }
  return iter->second;
}

