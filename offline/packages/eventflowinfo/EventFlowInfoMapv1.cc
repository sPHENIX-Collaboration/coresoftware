#include "EventFlowInfoMap.h"
#include "EventFlowInfoMapv1.h"


void EventFlowInfoMapv1::identify(std::ostream& os) const
{
  os << "EventFlowInfoMapv1: size = " << _map.size() << std::endl;
  for (const auto& m : _map)
  {
    m.second->identify(os);
  }
  return;
}

void EventFlowInfoMapv1::clear()
{
  for (const auto& iter : _map)
  {
    delete iter.second;
  }
  _map.clear();
  return;
}

const EventFlowInfo* EventFlowInfoMapv1::get(unsigned int id) const
{
  auto iter = _map.find(id);
  if (iter != _map.end())
  {
    return iter->second;
  }
  return nullptr;
}

EventFlowInfo* EventFlowInfoMapv1::get(unsigned int id)
{
  auto iter = _map.find(id);
  if (iter != _map.end())
  {
    return iter->second;
  }
  return nullptr;
}

const EventFlowInfo* EventFlowInfoMapv1::get(EventFlowInfo::EventFlowSrc type) const
{
  auto iter = _map.find(static_cast<unsigned int>(type));
  if (iter != _map.end())
  {
    return iter->second;
  }
  return nullptr;
}

EventFlowInfo* EventFlowInfoMapv1::get(EventFlowInfo::EventFlowSrc type)
{
  auto iter = _map.find(static_cast<unsigned int>(type));
  if (iter != _map.end())
  {
    return iter->second;
  }
  return nullptr;
}

EventFlowInfo* EventFlowInfoMapv1::insert(EventFlowInfo* clus, EventFlowInfo::EventFlowSrc type)
{
  unsigned int id = static_cast<unsigned int>(type);
  auto iter = _map.find(id);
  if (iter != _map.end())
  {
    delete iter->second; // delete existing object
  }
  _map[id] = clus; // insert new object
  return clus;
}

EventFlowInfoMap::ConstIter EventFlowInfoMapv1::find(EventFlowInfo::EventFlowSrc type) const
{
  return _map.find(static_cast<unsigned int>(type));
}

EventFlowInfoMap::Iter EventFlowInfoMapv1::find(EventFlowInfo::EventFlowSrc type)
{
  return _map.find(static_cast<unsigned int>(type));
}


