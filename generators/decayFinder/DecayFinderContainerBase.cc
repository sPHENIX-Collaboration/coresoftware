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

const DecayFinderContainerBase::Decay DecayFinderContainerBase::get(unsigned int) const
{
  return DummyDecay;
}

DecayFinderContainerBase::Decay DecayFinderContainerBase::get(unsigned int)
{
  return DummyDecay;
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
