/*****************/
/* Cameron Dean  */
/*   LANL 2021   */
/* cdean@bnl.gov */
/*****************/

#include "DecayFinderContainer_v1.h"

#include <map>       // for _Rb_tree_const_iterator, _Rb_tree_iterator
#include <ostream>   // for operator<<, endl, ostream, basic_ostream, bas...
#include <utility>   // for pair, make_pair
#include <vector>   // for vector

DecayFinderContainer_v1::DecayFinderContainer_v1()
{
}

DecayFinderContainer_v1::DecayFinderContainer_v1(const DecayFinderContainer_v1& decaymap):
DecayFinderContainerBase(decaymap)
{
  for (ConstIter iter = decaymap.begin(); iter != decaymap.end(); ++iter)
  {
    Decay decay = iter->second;
    m_decaymap.insert(std::make_pair(iter->first, decay));
  }
}

DecayFinderContainer_v1& DecayFinderContainer_v1::operator=(const DecayFinderContainer_v1& decaymap)
{
  Reset();
  for (ConstIter iter = decaymap.begin(); iter != decaymap.end(); ++iter)
  {
    Decay decay = iter->second;
    m_decaymap.insert(std::make_pair(iter->first, decay));
  }
  return *this;
}

DecayFinderContainer_v1::~DecayFinderContainer_v1()
{
  Reset();
}
