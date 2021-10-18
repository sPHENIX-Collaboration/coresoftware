/*****************/
/* Cameron Dean  */
/*   LANL 2021   */
/* cdean@bnl.gov */
/*****************/

#include "DecayFinderContainerBase.h"

#include <iterator>  // for reverse_iterator
#include <map>       // for _Rb_tree_const_iterator, _Rb_tree_iterator
#include <ostream>   // for operator<<, endl, ostream, basic_ostream, bas...
#include <utility>   // for pair, make_pair

DecayFinderContainerBase::DecayFinderContainerBase()
  : m_decaymap()
{
}

DecayFinderContainerBase::DecayFinderContainerBase(const DecayFinderContainerBase& decaymap)
  : m_decaymap()
{
  for (ConstIter iter = decaymap.begin(); iter != decaymap.end(); ++iter)
  {
    Decay decay = iter->second;
    m_decaymap.insert(std::make_pair(iter->first, decay));
  }
}

DecayFinderContainerBase& DecayFinderContainerBase::operator=(const DecayFinderContainerBase& decaymap)
{
  Reset();
  for (ConstIter iter = decaymap.begin(); iter != decaymap.end(); ++iter)
  {
    Decay decay = iter->second;
    m_decaymap.insert(std::make_pair(iter->first, decay));
  }
  return *this;
}

DecayFinderContainerBase::~DecayFinderContainerBase()
{
  Reset();
}

void DecayFinderContainerBase::Reset()
{
  m_decaymap.clear();
}

void DecayFinderContainerBase::identify(std::ostream& os) const
{
  os << "DecayFinderContainerBase: size = " << m_decaymap.size() << std::endl;
  return;
}

const DecayFinderContainerBase::Decay DecayFinderContainerBase::get(unsigned int id) const
{
  Decay dummyDecay = {{0, 0}};
  ConstIter iter = m_decaymap.find(id);
  if (iter == m_decaymap.end()) return dummyDecay;
  return iter->second;
}

DecayFinderContainerBase::Decay DecayFinderContainerBase::get(unsigned int id)
{
  Decay dummyDecay = {{0, 0}};
  Iter iter = m_decaymap.find(id);
  if (iter == m_decaymap.end()) return dummyDecay;
  return iter->second;
}

DecayFinderContainerBase::Decay DecayFinderContainerBase::insert(const Decay decay)
{
  unsigned int index = 0;
  if (!m_decaymap.empty()) index = m_decaymap.rbegin()->first + 1;
  m_decaymap.insert(std::make_pair(index, decay));
  return m_decaymap[index];
}

DecayFinderContainerBase::Map DecayFinderContainerBase::returnDecaysByPDGid(int PDGid)
{
  Map requiredDecays;

  for (Iter iter = m_decaymap.begin(); iter != m_decaymap.end(); ++iter)
    for (auto& [particleNumber, particleType] : iter->second)
      if (particleType == PDGid)
        requiredDecays.insert(std::make_pair(iter->first, iter->second));

  return requiredDecays;
}

size_t DecayFinderContainerBase::erase(unsigned int key)
{
  return m_decaymap.erase(key);
}
