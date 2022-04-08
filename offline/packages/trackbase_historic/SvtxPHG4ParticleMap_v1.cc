#include "SvtxPHG4ParticleMap_v1.h"

SvtxPHG4ParticleMap_v1::SvtxPHG4ParticleMap_v1() 
  : m_map()
{}

SvtxPHG4ParticleMap_v1::SvtxPHG4ParticleMap_v1(const SvtxPHG4ParticleMap_v1& map)
  : m_map()
{
  for (ConstIter iter = map.begin(); iter != map.end(); ++iter)
    {
      WeightedTruthTrackMap trackmap = iter->second;
      m_map.insert(std::make_pair(iter->first, trackmap));
    }
}


SvtxPHG4ParticleMap_v1& SvtxPHG4ParticleMap_v1::operator=(const SvtxPHG4ParticleMap_v1& map)
{
  clear();
  
  for(ConstIter iter = map.begin(); iter != map.end(); ++iter)
    {
      m_map.insert(std::make_pair(iter->first, iter->second));
    }

  return *this;
}

SvtxPHG4ParticleMap_v1::~SvtxPHG4ParticleMap_v1() 
{
}

void SvtxPHG4ParticleMap_v1::identify(std::ostream& os) const
{
  os << "SvtxPHG4ParticleMap_v1 size = " << m_map.size() << std::endl;
}

SvtxPHG4ParticleMap::WeightedTruthTrackMap SvtxPHG4ParticleMap_v1::insert(const unsigned int key, const SvtxPHG4ParticleMap::WeightedTruthTrackMap map)
{
  m_map.insert(std::make_pair(key, map)); 
  return map;
}
