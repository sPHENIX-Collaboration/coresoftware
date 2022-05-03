#include "PHG4ParticleSvtxMap_v1.h"

PHG4ParticleSvtxMap_v1::PHG4ParticleSvtxMap_v1()
  : m_map()
{
}

PHG4ParticleSvtxMap_v1::PHG4ParticleSvtxMap_v1(const PHG4ParticleSvtxMap_v1& map)
  : m_map()
{
  for (ConstIter iter = map.begin(); iter != map.end(); ++iter)
  {
    WeightedRecoTrackMap trackmap = iter->second;
    m_map.insert(std::make_pair(iter->first, trackmap));
  }
}

PHG4ParticleSvtxMap_v1& PHG4ParticleSvtxMap_v1::operator=(const PHG4ParticleSvtxMap_v1& map)
{
  Reset();

  for (ConstIter iter = map.begin(); iter != map.end(); ++iter)
  {
    m_map.insert(std::make_pair(iter->first, iter->second));
  }

  return *this;
}

PHG4ParticleSvtxMap_v1::~PHG4ParticleSvtxMap_v1()
{
  Reset();
}

void PHG4ParticleSvtxMap_v1::identify(std::ostream& os) const
{
  os << "PHG4ParticleSvtxMap_v1 size = " << m_map.size() << std::endl;

  for (const auto& [g4id, matchedRecoTracks] : m_map)
  {
    os << "G4Particle " << g4id << " has matched reco tracks " << std::endl;
    for (const auto& [weight, trackset] : matchedRecoTracks)
    {
      os << "  weight " << weight << " has reco tracks " << std::endl;
      for (const auto& track : trackset)
      {
        os << "    trackid : " << track << std::endl;
      }
    }
  }
}

PHG4ParticleSvtxMap::WeightedRecoTrackMap PHG4ParticleSvtxMap_v1::insert(const int key, const PHG4ParticleSvtxMap::WeightedRecoTrackMap map)
{
  m_map.insert(std::make_pair(key, map));
  return map;
}

const PHG4ParticleSvtxMap::WeightedRecoTrackMap & PHG4ParticleSvtxMap_v1::get(const int key) const
{
  const auto iter = m_map.find(key);
  if (iter == m_map.end())
  {
    return PHG4ParticleSvtxMap::get(key);
  }
  else
  {
    return iter->second;
  }
}
