/**
 * @file trackbase/CMFlashDifferenceContainerv1.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of CMFlashDifferenceContainerv1
 */
#include "CMFlashDifferenceContainerv1.h"
#include "CMFlashDifference.h"
#include "CMFlashDifferencev1.h"

#include <cstdlib>

void CMFlashDifferenceContainerv1::Reset()
{
  while (m_map.begin() != m_map.end())
  {
    delete m_map.begin()->second;
    m_map.erase(m_map.begin());
  }
  return;
}

void CMFlashDifferenceContainerv1::identify(std::ostream& os) const
{
  os << "-----CMFlashDifferenceContainerv1-----" << std::endl;
  ConstIterator iter;
  os << "Number of differences: " << size() << std::endl;
  for (iter = m_map.begin(); iter != m_map.end(); ++iter)
  {
    os << "key " << iter->first  << std::endl;
    (iter->second)->identify();
  }
  os << "------------------------------" << std::endl;
  return;
}

void
CMFlashDifferenceContainerv1::addDifferenceSpecifyKey(const unsigned int key, CMFlashDifference* newclus)
{
  auto ret = m_map.insert(std::make_pair(key, newclus));
  if ( !ret.second )
  {
    std::cout << "CMFlashDifferenceContainerv1::AddDifferenceSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
}

void CMFlashDifferenceContainerv1::removeDifference(unsigned int key)
{
  auto diff = findDifference(key);
  delete diff;
  m_map.erase(key); 
}

CMFlashDifferenceContainer::ConstRange
CMFlashDifferenceContainerv1::getDifferences() const
{ return std::make_pair(m_map.cbegin(), m_map.cend()); }

CMFlashDifference*
CMFlashDifferenceContainerv1::findDifference(unsigned int key) const
{
  auto it = m_map.find(key);
  return it == m_map.end() ? nullptr:it->second;
}

unsigned int CMFlashDifferenceContainerv1::size() const
{
  return m_map.size();
}
