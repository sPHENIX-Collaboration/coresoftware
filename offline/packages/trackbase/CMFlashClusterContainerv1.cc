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
  for( auto&& [key, cluster]:m_clusmap )
  { delete cluster; }
  
  m_clusmap.clear();
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

void CMFlashClusterContainerv1::addClusterSpecifyKey(const unsigned int key, CMFlashCluster* newclus)
{
  auto ret = m_clusmap.insert(std::make_pair(key, newclus));
  if ( !ret.second )
  {
    std::cout << "CMFlashClusterContainerv1::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
}

void CMFlashClusterContainerv1::removeCluster(unsigned int key)
{ 
  auto clus = findCluster(key);
  delete clus;

  m_clusmap.erase(key); 
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
