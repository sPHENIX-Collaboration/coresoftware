/**
 * @file trackbase/TrkrClusterContainerv1.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrClusterContainerv1
 */
#include "TrkrClusterContainerv1.h"
#include "TrkrCluster.h"
#include "TrkrClusterv2.h"

#include <cstdlib>

void TrkrClusterContainerv1::Reset()
{
  while (m_clusmap.begin() != m_clusmap.end())
  {
    delete m_clusmap.begin()->second;
    m_clusmap.erase(m_clusmap.begin());
  }
  return;
}

void TrkrClusterContainerv1::identify(std::ostream& os) const
{
  os << "-----TrkrClusterContainerv1-----" << std::endl;
  ConstIterator iter;
  os << "Number of clusters: " << size() << std::endl;
  for (iter = m_clusmap.begin(); iter != m_clusmap.end(); ++iter)
  {
    int layer = TrkrDefs::getLayer(iter->first);
    os << "clus key " << iter->first  << " layer " << layer << std::endl;
    (iter->second)->identify();
  }
  os << "------------------------------" << std::endl;
  return;
}

TrkrClusterContainerv1::ConstIterator
TrkrClusterContainerv1::addCluster(TrkrCluster* newclus)
{
  return addClusterSpecifyKey(newclus->getClusKey(), newclus);
}

TrkrClusterContainerv1::ConstIterator
TrkrClusterContainerv1::addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  auto ret = m_clusmap.insert(std::make_pair(key, newclus));
  if ( !ret.second )
  {
    std::cout << "TrkrClusterContainerv1::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
  else
  {
    return ret.first;
  }
}

void TrkrClusterContainerv1::removeCluster(TrkrDefs::cluskey key)
{ m_clusmap.erase(key); }

void TrkrClusterContainerv1::removeCluster(TrkrCluster *clus)
{ removeCluster( clus->getClusKey() ); }

TrkrClusterContainerv1::Iterator
TrkrClusterContainerv1::findOrAddCluster(TrkrDefs::cluskey key)
{
  auto it = m_clusmap.lower_bound( key );
  if (it == m_clusmap.end()|| (key < it->first ))
  {
    // add new cluster and set its key
    it = m_clusmap.insert(it, std::make_pair(key, new TrkrClusterv2()));
    it->second->setClusKey(key);
  }
  return it;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainerv1::getClusters() const
{ return std::make_pair(m_clusmap.cbegin(), m_clusmap.cend()); }

TrkrCluster*
TrkrClusterContainerv1::findCluster(TrkrDefs::cluskey key) const
{
  auto it = m_clusmap.find(key);
  return it == m_clusmap.end() ? nullptr:it->second;
}

unsigned int TrkrClusterContainerv1::size() const
{
  return m_clusmap.size();
}
