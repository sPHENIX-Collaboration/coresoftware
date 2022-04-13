#include "TrackSeedContainer_v1.h"

#include "TrackSeed.h"

#include <phool/PHObject.h>

#include <vector>
#include <iostream>

TrackSeedContainer_v1::TrackSeedContainer_v1() :
  m_seeds() {}

TrackSeedContainer_v1::TrackSeedContainer_v1(const TrackSeedContainer_v1& seeds) : m_seeds()
{
  for(const TrackSeed* seed : seeds)
    {
      TrackSeed* newseed = dynamic_cast<TrackSeed*>(seed->CloneMe());
      m_seeds.push_back(newseed);
    }
}

TrackSeedContainer_v1& TrackSeedContainer_v1::operator=(const TrackSeedContainer_v1& seedContainer)
{
  Reset();
  for(const TrackSeed* seed : seedContainer) 
    {
      TrackSeed* newseed = dynamic_cast<TrackSeed*>(seed->CloneMe());
      m_seeds.push_back(newseed);
    }

  return *this;
}

void TrackSeedContainer_v1::Reset() 
{
  for(const TrackSeed* seed : m_seeds)
    {
      delete seed;
    }

  m_seeds.clear();
}

// cppcheck-suppress virtualCallInConstructor
TrackSeedContainer_v1::~TrackSeedContainer_v1()
{
  Reset();
}

void TrackSeedContainer_v1::identify(std::ostream& os) const
{
  os << "TrackSeedContainer_v1 size is " << m_seeds.size()
     << std::endl; 
}

const TrackSeed* TrackSeedContainer_v1::get(const std::size_t key) const
{
  ConstIter iter = m_seeds.begin() + key;
  if(iter == m_seeds.end()) { return nullptr; }
  return *iter;
}

TrackSeed* TrackSeedContainer_v1::get(const std::size_t key) 
{
  Iter iter = m_seeds.begin() + key;
  if(iter == m_seeds.end()) { return nullptr; }
  return *iter;
}

TrackSeed* TrackSeedContainer_v1::insert(const TrackSeed* seed)
{
  m_seeds.push_back(dynamic_cast<TrackSeed*>(seed->CloneMe()));
  Iter iter = m_seeds.end() - 1;
  return *iter;
}
