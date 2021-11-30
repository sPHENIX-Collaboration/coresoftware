#include "GlobalVertex.h"

std::map<GlobalVertex::VTXTYPE, unsigned int> DummyGlobalVertex;

GlobalVertex::ConstVtxIter GlobalVertex::begin_vtxids() const
{
  return DummyGlobalVertex.end();
}

GlobalVertex::ConstVtxIter GlobalVertex::find_vtxids(VTXTYPE /*type*/) const
{
  return DummyGlobalVertex.end();
}

GlobalVertex::ConstVtxIter GlobalVertex::end_vtxids() const
{
  return DummyGlobalVertex.end();
}

GlobalVertex::VtxIter GlobalVertex::begin_vtxids()
{
  return DummyGlobalVertex.end();
}

GlobalVertex::VtxIter GlobalVertex::find_vtxids(VTXTYPE /*type*/)
{
  return DummyGlobalVertex.end();
}

GlobalVertex::VtxIter GlobalVertex::end_vtxids()
{
  return DummyGlobalVertex.end();
}
