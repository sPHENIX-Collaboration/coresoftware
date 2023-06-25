#include "Jet.h"

#include <iostream>

Jet::typ_comp_ids DummyJetMap;

void Jet::identify(std::ostream& os) const
{
  os << "---Jet-----------------------" << std::endl;
  return;
}

Jet::ConstIter Jet::begin_comp() const
{
  return DummyJetMap.end();
}

Jet::ConstIter Jet::lower_bound_comp(Jet::SRC /*source*/) const
{
  return DummyJetMap.end();
}

Jet::ConstIter Jet::upper_bound_comp(Jet::SRC /*source*/) const
{
  return DummyJetMap.end();
}

Jet::ConstIter Jet::find(Jet::SRC /*source*/) const
{
  return DummyJetMap.end();
}

Jet::ConstIter Jet::end_comp() const
{
  return DummyJetMap.end();
}

Jet::Iter Jet::begin_comp()
{
  return DummyJetMap.end();
}

Jet::Iter Jet::lower_bound_comp(Jet::SRC /*source*/)
{
  return DummyJetMap.end();
}

Jet::Iter Jet::upper_bound_comp(Jet::SRC /*source*/)
{
  return DummyJetMap.end();
}

Jet::Iter Jet::find(Jet::SRC /*source*/)
{
  return DummyJetMap.end();
}

Jet::Iter Jet::end_comp()
{
  return DummyJetMap.end();
}
