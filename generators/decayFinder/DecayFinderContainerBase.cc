/*****************/
/* Cameron Dean  */
/*   LANL 2021   */
/* cdean@bnl.gov */
/*****************/

#include "DecayFinderContainerBase.h"

//#include <iterator>  // for reverse_iterator
#include <map>       // for _Rb_tree_const_iterator, _Rb_tree_iterator
#include <ostream>   // for operator<<, endl, ostream, basic_ostream, bas...
#include <utility>   // for pair, make_pair

DecayFinderContainerBase::Decay DummyDecay;
DecayFinderContainerBase::Map DummyMap;

void DecayFinderContainerBase::Reset()
{
  DummyMap.clear();
}

void DecayFinderContainerBase::clear()
{
  DummyMap.clear();
}

const DecayFinderContainerBase::Decay DecayFinderContainerBase::get(unsigned int) const
{
  return DummyDecay;
}

DecayFinderContainerBase::Decay DecayFinderContainerBase::get(unsigned int)
{
  return DummyDecay;
}

DecayFinderContainerBase::ConstIter DecayFinderContainerBase::begin() const
{
  return DummyMap.end();
}

DecayFinderContainerBase::ConstIter DecayFinderContainerBase::find(unsigned int) const
{
  return DummyMap.end();
}

DecayFinderContainerBase::ConstIter DecayFinderContainerBase::end() const
{
  return DummyMap.end();
}


DecayFinderContainerBase::Iter DecayFinderContainerBase::begin()
{
  return DummyMap.end();
}

DecayFinderContainerBase::Iter DecayFinderContainerBase::find(unsigned int)
{
  return DummyMap.end();
}

DecayFinderContainerBase::Iter DecayFinderContainerBase::end()
{
  return DummyMap.end();
}

DecayFinderContainerBase::Decay DecayFinderContainerBase::insert(const Decay)
{
  return DummyDecay;
}

DecayFinderContainerBase::Map DecayFinderContainerBase::returnDecaysByPDGid(int)
{
  return DummyMap;
}

size_t DecayFinderContainerBase::erase(unsigned int key)
{
  return DummyMap.erase(key);
}
