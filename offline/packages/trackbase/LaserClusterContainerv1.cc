/**
 * @file trackbase/LaserClusterContainerv1.cc
 * @author Ben Kimelman
 * @date February 2024
 * @brief Implementation of LaserClusterContainerv1
 */
#include "LaserClusterContainerv1.h"
#include "LaserCluster.h"

#include <cstdlib>

void LaserClusterContainerv1::Reset()
{
  for( auto&& [key, cluster]:m_clusmap )
  { delete cluster; }
  
  m_clusmap.clear();
}

void LaserClusterContainerv1::identify(std::ostream& os) const
{
  os << "-----LaserClusterContainerv1-----" << std::endl;
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

void LaserClusterContainerv1::addClusterSpecifyKey(const unsigned int key, LaserCluster* newclus)
{
  auto ret = m_clusmap.insert(std::make_pair(key, newclus));
  if ( !ret.second )
  {
    std::cout << "LaserClusterContainerv1::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
}

void LaserClusterContainerv1::removeCluster(unsigned int key)
{ 
  auto clus = findCluster(key);
  delete clus;

  m_clusmap.erase(key); 
}

LaserClusterContainer::ConstRange
LaserClusterContainerv1::getClusters() const
{ return std::make_pair(m_clusmap.cbegin(), m_clusmap.cend()); }

LaserCluster*
LaserClusterContainerv1::findCluster(unsigned int key) const
{
  auto it = m_clusmap.find(key);
  return it == m_clusmap.end() ? nullptr:it->second;
}

unsigned int LaserClusterContainerv1::size() const
{
  return m_clusmap.size();
}
