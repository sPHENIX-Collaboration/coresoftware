/**
 * @file trackbase/TrkrClusterContainer.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrClusterContainer
 */
#include "TrkrClusterContainer.h"
#include "TrkrClusterv1.h"

#include <cstdlib>

TrkrClusterContainer::TrkrClusterContainer()
{
}

void TrkrClusterContainer::Reset()
{
  while (m_clusmap.begin() != m_clusmap.end())
  {
    delete m_clusmap.begin()->second;
    m_clusmap.erase(m_clusmap.begin());
  }
  return;
}

void TrkrClusterContainer::identify(std::ostream& os) const
{
  os << "-----TrkrClusterContainer-----" << std::endl;
  ConstIterator iter;
  os << "Number of clusters: " << size() << std::endl;
  for (iter = m_clusmap.begin(); iter != m_clusmap.end(); ++iter)
  {
    os << "clus key 0x" << std::hex << iter->first << std::dec << std::endl;
    (iter->second)->identify();
  }
  os << "------------------------------" << std::endl;
  return;
}

TrkrClusterContainer::ConstIterator
TrkrClusterContainer::addCluster(TrkrCluster* newclus)
{
  return addClusterSpecifyKey(newclus->getClusKey(), newclus);
}

TrkrClusterContainer::ConstIterator
TrkrClusterContainer::addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  auto ret = m_clusmap.insert(std::make_pair(key, newclus));
  if ( !ret.second )
  {
    std::cout << "TrkrClusterContainer::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
  else
  {
    return ret.first;
  }
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::getClusters(const TrkrDefs::TrkrId trackerid) const
{
  // TrkrDefs::cluskey tmp = trackerid;
  // TrkrDefs::cluskey keylow = tmp << TrackerDefs::bitshift_trackerid;
  // TrkrDefs::cluskey keyup = ((tmp + 1)<< TrackerDefs::bitshift_trackerid) -1 ;
  //   std::cout << "keylow: 0x" << hex << keylow << dec << std::endl;
  //   std::cout << "keyup: 0x" << hex << keyup << dec << std::endl;

  TrkrDefs::cluskey keylo = TrkrDefs::getClusKeyLo(trackerid);
  TrkrDefs::cluskey keyhi = TrkrDefs::getClusKeyHi(trackerid);

  ConstRange retpair;
  retpair.first = m_clusmap.lower_bound(keylo);
  retpair.second = m_clusmap.upper_bound(keyhi);
  return retpair;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::getClusters(const TrkrDefs::TrkrId trackerid, const char layer) const
{
  TrkrDefs::cluskey keylo = TrkrDefs::getClusKeyLo(trackerid, layer);
  TrkrDefs::cluskey keyhi = TrkrDefs::getClusKeyHi(trackerid, layer);

  ConstRange retpair;
  retpair.first = m_clusmap.lower_bound(keylo);
  retpair.second = m_clusmap.upper_bound(keyhi);
  return retpair;
}

TrkrClusterContainer::ConstRange
TrkrClusterContainer::getClusters(void) const
{
  return std::make_pair(m_clusmap.begin(), m_clusmap.end());
}

TrkrClusterContainer::Iterator
TrkrClusterContainer::findOrAddCluster(TrkrDefs::cluskey key)
{
  TrkrClusterContainer::Iterator it = m_clusmap.find(key);
  if (it == m_clusmap.end())
  {
    // add new cluster and set its key
    auto ret = m_clusmap.insert(std::make_pair(key, new TrkrClusterv1()));
    (ret.first->second)->setClusKey(key);
    it = ret.first;
  }
  return it;
}

TrkrCluster*
TrkrClusterContainer::findCluster(TrkrDefs::cluskey key)
{
  TrkrClusterContainer::Iterator it = m_clusmap.find(key);

  if (it != m_clusmap.end())
  {
    return it->second;
  }

  return 0;
}
