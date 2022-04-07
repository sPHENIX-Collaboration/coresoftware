#include "PHG4ParticleSvtxMap.h"

PHG4ParticleSvtxMap::Map DummyPHG4ParticleSvtxMap;
PHG4ParticleSvtxMap::WeightedRecoTrackMap emptyRecoMap;

const PHG4ParticleSvtxMap::WeightedRecoTrackMap PHG4ParticleSvtxMap::get(unsigned int) const 
{
  return emptyRecoMap;
}

PHG4ParticleSvtxMap::WeightedRecoTrackMap PHG4ParticleSvtxMap::get(unsigned int) 
{
  return emptyRecoMap;
}

PHG4ParticleSvtxMap::ConstIter PHG4ParticleSvtxMap::insert(const WeightedRecoTrackMap)
{
  return DummyPHG4ParticleSvtxMap.end();
}

PHG4ParticleSvtxMap::ConstIter PHG4ParticleSvtxMap::begin() const
{
  return DummyPHG4ParticleSvtxMap.end();
}

PHG4ParticleSvtxMap::ConstIter PHG4ParticleSvtxMap::find(unsigned int /*idkey*/) const
{
  return DummyPHG4ParticleSvtxMap.end();
}

PHG4ParticleSvtxMap::ConstIter PHG4ParticleSvtxMap::end() const
{
  return DummyPHG4ParticleSvtxMap.end();
}


PHG4ParticleSvtxMap::Iter PHG4ParticleSvtxMap::begin()
{
  return DummyPHG4ParticleSvtxMap.end();
}

PHG4ParticleSvtxMap::Iter PHG4ParticleSvtxMap::find(unsigned int /*idkey*/)
{
  return DummyPHG4ParticleSvtxMap.end();
}

PHG4ParticleSvtxMap::Iter PHG4ParticleSvtxMap::end()
{
  return DummyPHG4ParticleSvtxMap.end();
}
