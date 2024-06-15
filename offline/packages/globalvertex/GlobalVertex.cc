#include "GlobalVertex.h"

std::map<GlobalVertex::VTXTYPE, unsigned int> DummyGlobalVertex;
std::map<GlobalVertex::VTXTYPE, GlobalVertex::VertexVector> dumVertexVector;

GlobalVertex::ConstVertexIter GlobalVertex::begin_vertexes() const
{
  return dumVertexVector.end();
}
GlobalVertex::ConstVertexIter GlobalVertex::find_vertexes(VTXTYPE /*unused*/) const
{
  return dumVertexVector.end();
}
GlobalVertex::ConstVertexIter GlobalVertex::end_vertexes() const
{
  return dumVertexVector.end();
}
GlobalVertex::VertexIter GlobalVertex::begin_vertexes()
{
  return dumVertexVector.end();
}
GlobalVertex::VertexIter GlobalVertex::find_vertexes(VTXTYPE /*unused*/)
{
  return dumVertexVector.end();
}
GlobalVertex::VertexIter GlobalVertex::end_vertexes()
{
  return dumVertexVector.end();
}
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
