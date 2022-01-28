/**
 * @file trackbase/CMFlashClusterContainerv1.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of CMFlashClusterContainerv1
 */
#include "CMFlashClusterContainerv1.h"
#include "CMFlashCluster.h"
#include "CMFlashClusterv1.h"

#include <cstdlib>

void CMFlashClusterContainerv1::Reset()
{
  while (m_clusmap.begin() != m_clusmap.end())
  {
    delete m_clusmap.begin()->second;
    m_clusmap.erase(m_clusmap.begin());
  }
  return;
}

void CMFlashClusterContainerv1::identify(std::ostream& os) const
{
  os << "-----CMFlashClusterContainerv1-----" << std::endl;
  ConstIterator iter;
  os << "Number of clusters: " << size() << std::endl;
  for (iter = m_clusmap.begin(); iter != m_clusmap.end(); ++iter)
  {
    os << "clus key " << iter->first  << std::endl;
    (iter->second)->identify();
  }
  os << "------------------------------" << std::endl;
  return;
}

CMFlashClusterContainerv1::ConstIterator
CMFlashClusterContainerv1::addCluster(CMFlashCluster* newclus)
{
  return addClusterSpecifyKey(newclus->getClusKey(), newclus);
}

CMFlashClusterContainerv1::ConstIterator
CMFlashClusterContainerv1::addClusterSpecifyKey(const unsigned int key, CMFlashCluster* newclus)
{
  auto ret = m_clusmap.insert(std::make_pair(key, newclus));
  if ( !ret.second )
  {
    std::cout << "CMFlashClusterContainerv1::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
  else
  {
    return ret.first;
  }
}

void CMFlashClusterContainerv1::removeCluster(unsigned int key)
{ m_clusmap.erase(key); }

void CMFlashClusterContainerv1::removeCluster(CMFlashCluster *clus)
{ removeCluster( clus->getClusKey() ); }

CMFlashClusterContainerv1::Iterator
CMFlashClusterContainerv1::findOrAddCluster(unsigned int key)
{
  auto it = m_clusmap.lower_bound( key );
  if (it == m_clusmap.end()|| (key < it->first ))
  {
    // add new cluster and set its key
    it = m_clusmap.insert(it, std::make_pair(key, new CMFlashClusterv1()));
    it->second->setClusKey(key);
  }
  return it;
}

CMFlashClusterContainer::ConstRange
CMFlashClusterContainerv1::getClusters() const
{ return std::make_pair(m_clusmap.cbegin(), m_clusmap.cend()); }

CMFlashCluster*
CMFlashClusterContainerv1::findCluster(unsigned int key) const
{
  auto it = m_clusmap.find(key);
  return it == m_clusmap.end() ? nullptr:it->second;
}

unsigned int CMFlashClusterContainerv1::size() const
{
  return m_clusmap.size();
}
