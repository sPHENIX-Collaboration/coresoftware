/**
 * @file trackbase/TrkrClusterContainerv3.cc
 * @author D. McGlinchey, Hugo Pereira Da Costa
 * @date June 2018
 * @brief Implementation of TrkrClusterContainerv3
 */
#include "TrkrClusterContainerv3.h"
#include "TrkrCluster.h"
#include "TrkrClusterv3.h"
#include "TrkrDefs.h"

#include <cstdlib>

namespace
{
  TrkrClusterContainer::Map dummy_map;
}

//_________________________________________________________________
void TrkrClusterContainerv3::Reset()
{
  // delete all clusters
  for( const auto& map_pair:m_clusmap )
    for( const auto& pair:map_pair.second )
  { delete pair.second; }

  // clear the maps
  m_clusmap.clear();
}

//_________________________________________________________________
void TrkrClusterContainerv3::identify(std::ostream& os) const
{
  os << "-----TrkrClusterContainerv3-----" << std::endl;
  os << "Number of clusters: " << size() << std::endl;

  for( const auto& map_pair:m_clusmap )
  {

    const unsigned int layer = TrkrDefs::getLayer(map_pair.first);
    std::cout << "layer: " << layer << " hitsetkey: " << map_pair.first << std::endl;

    for( const auto& pair:map_pair.second )
    {
      int layer = TrkrDefs::getLayer(pair.first);
      os << "clus key " << pair.first  << " layer " << layer << std::endl;
      (pair.second)->identify();
    }
  }

  os << "------------------------------" << std::endl;
}

//_________________________________________________________________
void TrkrClusterContainerv3::removeCluster(TrkrDefs::cluskey key)
{
  // get hitset key from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey( key );
  
  // find relevant cluster map if any and remove corresponding cluster
  auto iter = m_clusmap.find( hitsetkey );
  if( iter != m_clusmap.end() ) iter->second.erase( key );
}
  
//_________________________________________________________________
void TrkrClusterContainerv3::removeCluster(TrkrCluster *clus)
{ removeCluster( clus->getClusKey() ); }

//_________________________________________________________________
TrkrClusterContainerv3::ConstIterator
TrkrClusterContainerv3::addCluster(TrkrCluster* newclus)
{ return addClusterSpecifyKey(newclus->getClusKey(), newclus); }

//_________________________________________________________________
TrkrClusterContainerv3::ConstIterator
TrkrClusterContainerv3::addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey( key );

  // find relevant cluster map or create one if not found
  Map& map = m_clusmap[hitsetkey];
  const auto ret = map.insert(std::make_pair(key, newclus));
  if ( !ret.second )
  {
    std::cout << "TrkrClusterContainerv3::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  } else {
    ret.first->second->setClusKey( key );
    return ret.first;
  }
}

//_________________________________________________________________
TrkrClusterContainerv3::ConstRange
TrkrClusterContainerv3::getClusters(TrkrDefs::hitsetkey hitsetkey) const
{
  // find relevant association map
  const auto iter = m_clusmap.find(hitsetkey);
  if( iter != m_clusmap.end() ) 
  {
    return std::make_pair( iter->second.cbegin(), iter->second.cend() );
  } else { 
    return std::make_pair( dummy_map.cbegin(), dummy_map.cend() );
  }
}

//_________________________________________________________________
TrkrClusterContainerv3::Map*
TrkrClusterContainerv3::getClusterMap(TrkrDefs::hitsetkey hitsetkey)
{ return &m_clusmap[hitsetkey]; }
  
//_________________________________________________________________
TrkrClusterContainerv3::Iterator
TrkrClusterContainerv3::findOrAddCluster(TrkrDefs::cluskey key)
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey( key );

  // find relevant cluster map or create one if not found
  Map& map = m_clusmap[hitsetkey];
  auto it = map.lower_bound(key);
  if( it == map.end() || (key<it->first) )
  {
    // add new cluster and set its key
    it = map.insert(it, std::make_pair(key, new TrkrClusterv3()));
    it->second->setClusKey( key );
  }
  
  return it;
}

//_________________________________________________________________
TrkrCluster* TrkrClusterContainerv3::findCluster(TrkrDefs::cluskey key) const
{
  
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey( key );

  const auto map_iter = m_clusmap.find(hitsetkey);
  if( map_iter != m_clusmap.end() ) 
  {
    const auto clus_iter = map_iter->second.find( key );
    if( clus_iter !=  map_iter->second.end() )
    { 
      return clus_iter->second;
    } else {
      return nullptr;
    }
  } else {
    return nullptr;
  }
}

//_________________________________________________________________
unsigned int TrkrClusterContainerv3::size(void) const
{
  unsigned int size = 0;
  for( const auto& map_pair:m_clusmap )
  { size += map_pair.second.size(); }
  
  return size;
}
