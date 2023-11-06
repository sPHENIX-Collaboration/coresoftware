#include "EventplaneinfoMap.h"

std::map<unsigned int, Eventplaneinfo*> DummyEventplaneinfoMap;

EventplaneinfoMap::ConstIter EventplaneinfoMap::begin() const
{
  return DummyEventplaneinfoMap.end();
}

EventplaneinfoMap::ConstIter EventplaneinfoMap::find(unsigned int /*idkey*/) const
{
  return DummyEventplaneinfoMap.end();
}

EventplaneinfoMap::ConstIter EventplaneinfoMap::end() const
{
  return DummyEventplaneinfoMap.end();
}

EventplaneinfoMap::Iter EventplaneinfoMap::begin()
{
  return DummyEventplaneinfoMap.end();
}

EventplaneinfoMap::Iter EventplaneinfoMap::find(unsigned int /*idkey*/)
{
  return DummyEventplaneinfoMap.end();
}

EventplaneinfoMap::Iter EventplaneinfoMap::end()
{
  return DummyEventplaneinfoMap.end();
}
