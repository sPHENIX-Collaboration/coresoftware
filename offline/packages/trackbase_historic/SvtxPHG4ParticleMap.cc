#include "SvtxPHG4ParticleMap.h"

SvtxPHG4ParticleMap::Map DummySvtxPHG4ParticleMap;
SvtxPHG4ParticleMap::WeightedTruthTrackMap emptyTruthMap;

const SvtxPHG4ParticleMap::WeightedTruthTrackMap & SvtxPHG4ParticleMap::get(const unsigned int) const
{
  return emptyTruthMap;
}

SvtxPHG4ParticleMap::WeightedTruthTrackMap & SvtxPHG4ParticleMap::get(const unsigned int)
{
  return emptyTruthMap;
}

SvtxPHG4ParticleMap::WeightedTruthTrackMap SvtxPHG4ParticleMap::insert(const unsigned int, const WeightedTruthTrackMap)
{
  return emptyTruthMap;
}

SvtxPHG4ParticleMap::ConstIter SvtxPHG4ParticleMap::begin() const
{
  return DummySvtxPHG4ParticleMap.end();
}

SvtxPHG4ParticleMap::ConstIter SvtxPHG4ParticleMap::find(const unsigned int) const
{
  return DummySvtxPHG4ParticleMap.end();
}

SvtxPHG4ParticleMap::ConstIter SvtxPHG4ParticleMap::end() const
{
  return DummySvtxPHG4ParticleMap.end();
}

SvtxPHG4ParticleMap::Iter SvtxPHG4ParticleMap::begin()
{
  return DummySvtxPHG4ParticleMap.end();
}

SvtxPHG4ParticleMap::Iter SvtxPHG4ParticleMap::find(unsigned int)
{
  return DummySvtxPHG4ParticleMap.end();
}

SvtxPHG4ParticleMap::Iter SvtxPHG4ParticleMap::end()
{
  return DummySvtxPHG4ParticleMap.end();
}
