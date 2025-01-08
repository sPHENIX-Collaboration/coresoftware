#include "JetMap.h"

#include <ostream>  // for operator<<, endl, ostream, basic_ostream

JetMap::typ_JetMap DummyJetMapMap;
std::set<Jet::SRC> DummyJetSourceSet;

void JetMap::identify(std::ostream& os) const
{
  os << "JetMap" << std::endl;
  return;
}

JetMap::ConstSrcIter JetMap::begin_src() const
{
  return DummyJetSourceSet.end();
}

JetMap::ConstSrcIter JetMap::find_src(Jet::SRC /*src*/) const
{
  return DummyJetSourceSet.end();
}

JetMap::ConstSrcIter JetMap::end_src() const
{
  return DummyJetSourceSet.end();
}

JetMap::SrcIter JetMap::begin_src()
{
  return DummyJetSourceSet.end();
}

JetMap::SrcIter JetMap::find_src(Jet::SRC /*src*/)
{
  return DummyJetSourceSet.end();
}

JetMap::SrcIter JetMap::end_src()
{
  return DummyJetSourceSet.end();
}

JetMap::ConstIter JetMap::begin() const
{
  return DummyJetMapMap.end();
}

JetMap::ConstIter JetMap::find(unsigned int /*idkey*/) const
{
  return DummyJetMapMap.end();
}

JetMap::ConstIter JetMap::end() const
{
  return DummyJetMapMap.end();
}

JetMap::Iter JetMap::begin()
{
  return DummyJetMapMap.end();
}

JetMap::Iter JetMap::find(unsigned int /*idkey*/)
{
  return DummyJetMapMap.end();
}

JetMap::Iter JetMap::end()
{
  return DummyJetMapMap.end();
}
