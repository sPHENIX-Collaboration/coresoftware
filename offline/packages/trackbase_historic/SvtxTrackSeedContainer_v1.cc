#include "SvtxTrackSeedContainer_v1.h"

#include "SvtxTrackSeed.h"

#include <phool/PHObject.h>

#include <vector>
#include <iostream>

SvtxTrackSeedContainer_v1::SvtxTrackSeedContainer_v1() :
  m_seeds() {}

SvtxTrackSeedContainer_v1::SvtxTrackSeedContainer_v1(const SvtxTrackSeedContainer_v1& seeds) : m_seeds()
{
  for(const SvtxTrackSeed* seed : seeds)
    {
      SvtxTrackSeed* newseed = dynamic_cast<SvtxTrackSeed*>(seed->CloneMe());
      m_seeds.push_back(newseed);
    }
}

SvtxTrackSeedContainer_v1& SvtxTrackSeedContainer_v1::operator=(const SvtxTrackSeedContainer_v1& seedContainer)
{
  Reset();
  for(const SvtxTrackSeed* seed : seedContainer) 
    {
      SvtxTrackSeed* newseed = dynamic_cast<SvtxTrackSeed*>(seed->CloneMe());
      m_seeds.push_back(newseed);
    }

  return *this;
}

SvtxTrackSeedContainer_v1::~SvtxTrackSeedContainer_v1()
{
  Reset();
}

void SvtxTrackSeedContainer_v1::Reset() 
{
  for(const SvtxTrackSeed* seed : m_seeds)
    {
      delete seed;
    }

  m_seeds.clear();
}

void SvtxTrackSeedContainer_v1::identify(std::ostream& os) const
{
  os << "SvtxTrackSeedContainer_v1 size is " << m_seeds.size()
     << std::endl; 
}

const SvtxTrackSeed* SvtxTrackSeedContainer_v1::get(const std::size_t key) const
{
  ConstIter iter = m_seeds.begin() + key;
  if(iter == m_seeds.end()) { return nullptr; }
  return *iter;
}

SvtxTrackSeed* SvtxTrackSeedContainer_v1::get(const std::size_t key) 
{
  Iter iter = m_seeds.begin() + key;
  if(iter == m_seeds.end()) { return nullptr; }
  return *iter;
}

SvtxTrackSeed* SvtxTrackSeedContainer_v1::insert(const SvtxTrackSeed* seed)
{
  m_seeds.push_back(dynamic_cast<SvtxTrackSeed*>(seed->CloneMe()));
  Iter iter = m_seeds.end() - 1;
  return *iter;
}
