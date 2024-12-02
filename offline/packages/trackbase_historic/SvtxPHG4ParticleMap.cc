#include "SvtxPHG4ParticleMap.h"

SvtxPHG4ParticleMap::Map DummySvtxPHG4ParticleMap;
SvtxPHG4ParticleMap::WeightedTruthTrackMap emptyTruthMap;

const SvtxPHG4ParticleMap::WeightedTruthTrackMap& SvtxPHG4ParticleMap::get(const unsigned int /*unused*/) const
{
  return emptyTruthMap;
}

SvtxPHG4ParticleMap::WeightedTruthTrackMap& SvtxPHG4ParticleMap::get(const unsigned int /*unused*/)
{
  return emptyTruthMap;
}

SvtxPHG4ParticleMap::WeightedTruthTrackMap SvtxPHG4ParticleMap::insert(const unsigned int /*unused*/, const WeightedTruthTrackMap /*unused*/)
{
  return emptyTruthMap;
}

SvtxPHG4ParticleMap::ConstIter SvtxPHG4ParticleMap::begin() const
{
  return DummySvtxPHG4ParticleMap.end();
}

SvtxPHG4ParticleMap::ConstIter SvtxPHG4ParticleMap::find(const unsigned int /*unused*/) const
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

SvtxPHG4ParticleMap::Iter SvtxPHG4ParticleMap::find(unsigned int /*unused*/)
{
  return DummySvtxPHG4ParticleMap.end();
}

SvtxPHG4ParticleMap::Iter SvtxPHG4ParticleMap::end()
{
  return DummySvtxPHG4ParticleMap.end();
}
