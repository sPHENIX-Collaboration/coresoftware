#include "SvtxAlignmentStateMap_v1.h"

#include "SvtxAlignmentState.h"

#include <phool/PHObject.h>

SvtxAlignmentStateMap_v1::SvtxAlignmentStateMap_v1()
  : m_map()
{}

SvtxAlignmentStateMap_v1::SvtxAlignmentStateMap_v1(const SvtxAlignmentStateMap_v1& map)
  : m_map()
{
  for(auto [key, state] : map)
    {
      auto statecl = static_cast<SvtxAlignmentState*>(state->CloneMe());
      m_map.insert(std::make_pair(key, statecl));
    }
}

SvtxAlignmentStateMap_v1& SvtxAlignmentStateMap_v1::operator=(const SvtxAlignmentStateMap_v1& map)
{
  if (&map == this) 
    {
      return *this;
    }

  Reset();
  for(auto [key, state] : map)
    {
      auto clstate = static_cast<SvtxAlignmentState*>(state->CloneMe());
      m_map.insert(std::make_pair(key, clstate));
    }

  return *this;
}

SvtxAlignmentStateMap_v1::~SvtxAlignmentStateMap_v1()
{
  Reset();
}

void SvtxAlignmentStateMap_v1::Reset()
{
  for(auto [key, state] : m_map)
    {
      delete state;
    }
 
  m_map.clear();

}


void SvtxAlignmentStateMap_v1::identify(std::ostream& os) const
{
  os << "SvtxAlignmentStateMap_v1: size = " << m_map.size() << std::endl;
}

const SvtxAlignmentState* SvtxAlignmentStateMap_v1::get(TrkrDefs::cluskey ckey) const
{
  auto iter = m_map.find(ckey);
  if(iter == m_map.end()) 
    {
      return nullptr;
    }

  return iter->second;
}

SvtxAlignmentState* SvtxAlignmentStateMap_v1::get(TrkrDefs::cluskey ckey)
{
  auto iter = m_map.find(ckey);
  if(iter == m_map.end()) 
    {
      return nullptr;
    }

  return iter->second;
}

SvtxAlignmentState* SvtxAlignmentStateMap_v1::insertWithKey(TrkrDefs::cluskey ckey, SvtxAlignmentState* state)
{
  auto copy = static_cast<SvtxAlignmentState*> ( state->CloneMe());
  const auto result = m_map.insert(std::make_pair(ckey, copy));
  
  if(!result.second)
    {
      std::cout <<  "SvtxAlignmentStateMap_v1::insertWithKey - duplicated key, state not entered. " << std::endl;
      delete copy;
      return nullptr;
    }
  else
    {
      return copy;
    }
}
