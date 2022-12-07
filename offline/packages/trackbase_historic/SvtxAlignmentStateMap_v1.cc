#include "SvtxAlignmentStateMap_v1.h"

#include "SvtxAlignmentState.h"

#include <phool/PHObject.h>

SvtxAlignmentStateMap_v1::SvtxAlignmentStateMap_v1()
  : m_map()
{
}

SvtxAlignmentStateMap_v1::~SvtxAlignmentStateMap_v1()
{
  Reset();
}

void SvtxAlignmentStateMap_v1::Reset()
{
  for (auto [key, statevec] : m_map)
  {
    for (auto state : statevec)
    {
      delete state;
    }
    statevec.clear();
  }

  m_map.clear();
}

std::size_t SvtxAlignmentStateMap_v1::erase(unsigned int track)
{
  auto iter = m_map.find(track);
  if (iter != m_map.end())
  {
    auto statevec = iter->second;
    for (auto state : statevec)
    {
      delete state;
    }
  }

  return m_map.erase(track);
}

void SvtxAlignmentStateMap_v1::identify(std::ostream& os) const
{
  os << "SvtxAlignmentStateMap_v1: size = " << m_map.size() << std::endl;
}

const SvtxAlignmentStateMap::StateVec SvtxAlignmentStateMap_v1::get(unsigned int track) const
{
  auto iter = m_map.find(track);
  if (iter == m_map.end())
  {
    return SvtxAlignmentStateMap::StateVec();
  }

  return iter->second;
}

SvtxAlignmentStateMap::StateVec SvtxAlignmentStateMap_v1::get(unsigned int track)
{
  auto iter = m_map.find(track);
  if (iter == m_map.end())
  {
    return SvtxAlignmentStateMap::StateVec();
  }

  return iter->second;
}

SvtxAlignmentStateMap::StateVec SvtxAlignmentStateMap_v1::insertWithKey(unsigned int track, SvtxAlignmentStateMap::StateVec states)
{
  const auto result = m_map.insert(std::make_pair(track, states));

  if (!result.second)
  {
    std::cout << "SvtxAlignmentStateMap_v1::insertWithKey - duplicated key " << track << ", state not entered. " << std::endl;

    return SvtxAlignmentStateMap::StateVec();
  }
  else
  {
    return states;
  }
}
