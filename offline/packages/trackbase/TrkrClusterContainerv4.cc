/**
 * @file trackbase/TrkrClusterContainerv4.cc
 * @author D. McGlinchey, Hugo Pereira Da Costa
 * @date June 2018
 * @brief Implementation of TrkrClusterContainerv4
 */
#include "TrkrClusterContainerv4.h"
#include "TrkrCluster.h"
#include "TrkrClusterv3.h"
#include "TrkrDefs.h"

#include <cstdlib>

namespace
{
  TrkrClusterContainer::Map dummy_map;
}

//_________________________________________________________________
void TrkrClusterContainerv4::Reset()
{
  // delete all clusters
  for( auto&& [key,clus_vector]:m_clusmap )
  {
    for( auto&& cluster:clus_vector )
    { delete cluster; }
  }

  // clear the maps
  /* using swap ensures that the memory is properly de-allocated */
  {
    std::map<TrkrDefs::hitsetkey, Vector> empty;
    m_clusmap.swap( empty );
  }
  
  // also clear temporary map
  {
    Map empty;
    m_tmpmap.swap( empty );
  }
}


//_________________________________________________________________
void TrkrClusterContainerv4::identify(std::ostream& os) const
{
  os << "-----TrkrClusterContainerv4-----" << std::endl;
  os << "Number of clusters: " << size() << std::endl;

  for( const auto& [hitsetkey,clus_vector]:m_clusmap )
  {
    const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    std::cout << "layer: " << layer << " hitsetkey: " << hitsetkey << std::endl;

    for( const auto& cluster:clus_vector )
    { cluster->identify(); }
  }

  os << "------------------------------" << std::endl;
}

//_________________________________________________________________
void TrkrClusterContainerv4::removeCluster(TrkrDefs::cluskey key)
{
  // get hitset key from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey( key );
  
  // find relevant cluster map if any and remove corresponding cluster
  auto iter = m_clusmap.find( hitsetkey );
  if( iter != m_clusmap.end() )
  {
    // local reference to the vector
    auto& clus_vector = iter->second;
    
    // cluster index in vector
    const auto index = TrkrDefs::getClusIndex(key);

    // compare to vector size
    if( index < clus_vector.size() )
    {
      // delete corresponding element and set to null
      delete clus_vector[index];
      clus_vector[index] = nullptr;
    }
  }
}
  
//_________________________________________________________________
void TrkrClusterContainerv4::removeCluster(TrkrCluster *clus)
{ removeCluster( clus->getClusKey() ); }

//_________________________________________________________________
void
TrkrClusterContainerv4::addCluster(TrkrCluster* newclus)
{ return addClusterSpecifyKey(newclus->getClusKey(), newclus); }

//_________________________________________________________________
void
TrkrClusterContainerv4::addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey( key );

  // find relevant vector or create one if not found
  auto& clus_vector = m_clusmap[hitsetkey];
  
  // get cluster index in vector
  const auto index = TrkrDefs::getClusIndex( key );
  
  // compare index to vector size
  if( index < clus_vector.size() )
  {
    /* 
     * if index is already contained in vector, check corresponding element
     * and assign newclus if null
     * print error message and exit otherwise
     */
    if( !clus_vector[index] ) clus_vector[index] = newclus;
    else {
      std::cout << "TrkrClusterContainerv4::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
      exit(1);
    }
  } else if( index == clus_vector.size() ) {
    // if index matches the vector size, just push back the new cluster
    clus_vector.push_back( newclus );
  } else {
    // if index exceeds the vector size, resize cluster to the right size with nullptr, and assign
    clus_vector.resize( index+1, nullptr );
    clus_vector[index] = newclus;
  }
}

//_________________________________________________________________
TrkrClusterContainerv4::ConstRange
TrkrClusterContainerv4::getClusters(TrkrDefs::hitsetkey hitsetkey)
{
  // clear temporary map
  {
    Map empty;
    m_tmpmap.swap( empty );
  }

  // find relevant vector
  const auto iter = m_clusmap.find(hitsetkey);
  if( iter != m_clusmap.end() ) 
  {
    // copy content in temporary map
    for( const auto& cluster:iter->second )
    { if( cluster ) m_tmpmap.insert( m_tmpmap.end(), std::make_pair( cluster->getClusKey, cluster ) ); }
  }

  // return temporary map range
  return std::make_pair( m_tmpmap.cbegin(), m_tmpmap.cend() );
}

//_________________________________________________________________
TrkrClusterContainerv4::Vector*
TrkrClusterContainerv4::getClusterVector(TrkrDefs::hitsetkey hitsetkey)
{ return &m_clusmap[hitsetkey]; }
  
//_________________________________________________________________
TrkrCluster* TrkrClusterContainerv4::findCluster(TrkrDefs::cluskey key) const
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey( key );

  const auto map_iter = m_clusmap.find(hitsetkey);
  if( map_iter != m_clusmap.end() ) 
  {
    // local reference to vector
    const auto& clus_vector = map_iter->second;
    
    // get cluster position in vector
    const int index = TrkrDefs::getClusIndex( key );
    
    // compare to vector size
    if( index < clus_vector.size() ) return clus_vector[index];
    else return nullptr;
  } else return nullptr;
}

//_________________________________________________________________
unsigned int TrkrClusterContainerv4::size(void) const
{
  unsigned int size = 0;
  for( const auto& [hitsetkey,clus_vector]:m_clusmap )
  { size += clus_vector.size(); }
  return size;
}
