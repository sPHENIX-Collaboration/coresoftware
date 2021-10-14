/*****************/
/* Cameron Dean  */
/*   LANL 2021   */
/* cdean@bnl.gov */
/*****************/

#include "DecayFinderContainer.h"

#include <phool/PHObject.h>  // for PHObject

#include <iterator>  // for reverse_iterator
#include <map>       // for _Rb_tree_const_iterator, _Rb_tree_iterator
#include <ostream>   // for operator<<, endl, ostream, basic_ostream, bas...
#include <utility>   // for pair, make_pair

DecayFinderContainer::DecayFinderContainer()
  : m_decaymap()
{
}

DecayFinderContainer::DecayFinderContainer(const DecayFinderContainer& decaymap)
  : m_decaymap()
{
  for (ConstIter iter = decaymap.begin(); iter != decaymap.end(); ++iter)
  {
    Decay decay = iter->second;
    m_decaymap.insert(std::make_pair(iter->first, decay));
  }
}

DecayFinderContainer& DecayFinderContainer::operator=(const DecayFinderContainer& decaymap)
{
  Reset();
  for (ConstIter iter = decaymap.begin(); iter != decaymap.end(); ++iter)
  {
    Decay decay = iter->second;
    m_decaymap.insert(std::make_pair(iter->first, decay));
  }
  return *this;
}

DecayFinderContainer::~DecayFinderContainer()
{
  Reset();
}

void DecayFinderContainer::Reset()
{
  m_decaymap.clear();
}

void DecayFinderContainer::identify(std::ostream& os) const
{
  os << "DecayFinderContainer: size = " << m_decaymap.size() << std::endl;
  return;
}

const DecayFinderContainer::Decay DecayFinderContainer::get(unsigned int id) const
{
  Decay dummyDecay = {{0, 0}};
  ConstIter iter = m_decaymap.find(id);
  if (iter == m_decaymap.end()) return dummyDecay;
  return iter->second;
}

DecayFinderContainer::Decay DecayFinderContainer::get(unsigned int id)
{
  Decay dummyDecay = {{0, 0}};
  Iter iter = m_decaymap.find(id);
  if (iter == m_decaymap.end()) return dummyDecay;
  return iter->second;
}

DecayFinderContainer::Decay DecayFinderContainer::insert(const Decay decay)
{
  unsigned int index = 0;
  if (!m_decaymap.empty()) index = m_decaymap.rbegin()->first + 1;
  m_decaymap.insert(std::make_pair(index, decay));
  return m_decaymap[index];
}

DecayFinderContainer::Map DecayFinderContainer::returnDecaysByPDGid(int PDGid)
{
  Map requiredDecays;

  for (Iter iter = m_decaymap.begin(); iter != m_decaymap.end(); ++iter)
    for (auto& [particleNumber, particleType] : iter->second)
      if (particleType == PDGid)
        requiredDecays.insert(std::make_pair(iter->first, iter->second));

  return requiredDecays;
}

size_t DecayFinderContainer::erase(unsigned int key)
{
  return m_decaymap.erase(key);
}
