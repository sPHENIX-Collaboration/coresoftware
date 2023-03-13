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
  : m_decaymap()
{
}

DecayFinderContainer_v1::DecayFinderContainer_v1(const DecayFinderContainer_v1& decaymap)
  : m_decaymap()
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

void DecayFinderContainer_v1::Reset()
{
  m_decaymap.clear();
}

void DecayFinderContainer_v1::identify(std::ostream& os) const
{
  os << "DecayFinderContainer_v1: size = " << m_decaymap.size() << std::endl;
  return;
}

const DecayFinderContainerBase::Decay DecayFinderContainer_v1::get(unsigned int id) const
{
  std::pair<int, int> dummyPair = {0, 0};
  Decay dummyDecay = {{dummyPair, 0}};
  ConstIter iter = m_decaymap.find(id);
  if (iter == m_decaymap.end()) return dummyDecay;
  return iter->second;
}

DecayFinderContainerBase::Decay DecayFinderContainer_v1::get(unsigned int id)
{
  std::pair<int, int> dummyPair = {0, 0};
  Decay dummyDecay = {{dummyPair, 0}};
  Iter iter = m_decaymap.find(id);
  if (iter == m_decaymap.end()) return dummyDecay;
  return iter->second;
}

DecayFinderContainerBase::Decay DecayFinderContainer_v1::insert(const Decay decay)
{
  unsigned int index = 0;
  if (!m_decaymap.empty()) index = m_decaymap.rbegin()->first + 1;
  m_decaymap.insert(std::make_pair(index, decay));
  return m_decaymap[index];
}

DecayFinderContainerBase::Map DecayFinderContainer_v1::returnDecaysByPDGid(int PDGid)
{
  Map requiredDecays;

  for (Iter iter = m_decaymap.begin(); iter != m_decaymap.end(); ++iter)
    for (auto& [particleNumber, particleType] : iter->second)
      if (particleType == PDGid)
        requiredDecays.insert(std::make_pair(iter->first, iter->second));

  return requiredDecays;
}

size_t DecayFinderContainer_v1::erase(unsigned int key)
{
  return m_decaymap.erase(key);
}
