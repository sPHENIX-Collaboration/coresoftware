#include "PHG4Shower.h"

PHG4Shower::ParticleIdSet DummyParticleIdSet;
PHG4Shower::VertexIdSet DummyVertexIdSet;
PHG4Shower::HitIdMap DummyHitIdMap;

PHG4Shower::ParticleIdConstIter PHG4Shower::begin_g4particle_id() const
{
  return DummyParticleIdSet.end();
}

PHG4Shower::ParticleIdConstIter PHG4Shower::end_g4particle_id() const
{
  return DummyParticleIdSet.end();
}

PHG4Shower::ParticleIdIter PHG4Shower::begin_g4particle_id()
{
  return DummyParticleIdSet.end();
}

PHG4Shower::ParticleIdIter PHG4Shower::end_g4particle_id()
{
  return DummyParticleIdSet.end();
}

PHG4Shower::VertexIdConstIter PHG4Shower::begin_g4vertex_id() const
{
  return DummyVertexIdSet.end();
}

PHG4Shower::VertexIdConstIter PHG4Shower::end_g4vertex_id() const
{
  return DummyVertexIdSet.end();
}

PHG4Shower::VertexIdIter PHG4Shower::begin_g4vertex_id()
{
  return DummyVertexIdSet.end();
}

PHG4Shower::VertexIdIter PHG4Shower::end_g4vertex_id()
{
  return DummyVertexIdSet.end();
}

PHG4Shower::HitIdConstIter PHG4Shower::begin_g4hit_id() const
{
  return DummyHitIdMap.end();
}

PHG4Shower::HitIdConstIter PHG4Shower::find_g4hit_id(int /*volume*/) const
{
  return DummyHitIdMap.end();
}

PHG4Shower::HitIdConstIter PHG4Shower::end_g4hit_id() const
{
  return DummyHitIdMap.end();
}

PHG4Shower::HitIdIter PHG4Shower::begin_g4hit_id()
{
  return DummyHitIdMap.end();
}

PHG4Shower::HitIdIter PHG4Shower::find_g4hit_id(int /*volume*/)
{
  return DummyHitIdMap.end();
}

PHG4Shower::HitIdIter PHG4Shower::end_g4hit_id()
{
  return DummyHitIdMap.end();
}
