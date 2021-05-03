/**
 * @file trackbase/TrkrClusterContainerv1.cc
 * @author D. McGlinchey, Hugo Pereira Da Costa
 * @date June 2018
 * @brief Implementation of TrkrClusterContainerv1
 */
#include "TrkrClusterContainerv1.h"
#include "TrkrCluster.h"
#include "TrkrClusterv1.h"
#include "TrkrDefs.h"

#include <cstdlib>

namespace
{
  TrkrClusterContainer::Map dummy_map;
}

//_________________________________________________________________
void TrkrClusterContainerv1::Reset()
{
  for(unsigned int layer = 0;layer < max_layer; layer++)
    for(unsigned int phi_segment = 0;phi_segment < max_phisegment;phi_segment++)
    for(unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++)
  {
    for( const auto& pair:m_clusmap[layer][phi_segment][z_segment] )
    { delete pair.second; }

    m_clusmap[layer][phi_segment][z_segment].clear();
  }
}

//_________________________________________________________________
void TrkrClusterContainerv1::identify(std::ostream& os) const
{
  os << "-----TrkrClusterContainerv1-----" << std::endl;
  os << "Number of clusters: " << size() << std::endl;

  for(unsigned int layer = 0;layer < max_layer; layer++)
    for(unsigned int phi_segment = 0;phi_segment < max_phisegment;phi_segment++)
    for(unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++)
  {
    const auto iter = m_clusmap[layer][phi_segment][z_segment].begin();
    const int flayer = TrkrDefs::getLayer(iter->first);
    const unsigned int fsector = TrkrDefs::getPhiElement(iter->first);
    const unsigned int fside = TrkrDefs::getZElement(iter->first);
    std::cout << "layer: " << layer << " | " << flayer
      << " phi_seg: " << phi_segment << " | " << fsector
      << " z_seg: " << z_segment << " | " << fside
      << " nclu: " << m_clusmap[layer][phi_segment][z_segment].size()
      << std::endl;

    for( const auto& pair:m_clusmap[layer][phi_segment][z_segment] )
    {
      int layer = TrkrDefs::getLayer(pair.first);
      os << "clus key " << pair.first  << " layer " << layer << std::endl;
      (pair.second)->identify();
    }
  }

  os << "------------------------------" << std::endl;
}

//_________________________________________________________________
void TrkrClusterContainerv1::removeCluster(TrkrDefs::cluskey key)
{
  unsigned int layer = TrkrDefs::getLayer(key);
  unsigned int sector = TrkrDefs::getPhiElement(key);
  unsigned int side = TrkrDefs::getZElement(key);

  // bound check
  if( layer < max_layer && sector < max_phisegment && side < max_zsegment )
  {

    m_clusmap[layer][sector][side].erase(key);

  } else {

    std::cout
      << "TrkrClusterContainerv1::removeCluster - out of range access."
      << " layer: " << layer
      << " sector: " << sector
      << " side: " << side
      << std::endl;

  }
}

//_________________________________________________________________
void TrkrClusterContainerv1::removeCluster(TrkrCluster *clus)
{ removeCluster( clus->getClusKey() ); }

//_________________________________________________________________
TrkrClusterContainerv1::ConstIterator
TrkrClusterContainerv1::addCluster(TrkrCluster* newclus)
{ return addClusterSpecifyKey(newclus->getClusKey(), newclus); }

//_________________________________________________________________
TrkrClusterContainerv1::ConstIterator
TrkrClusterContainerv1::addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  unsigned int layer = TrkrDefs::getLayer(key);
  unsigned int sector = TrkrDefs::getPhiElement(key);
  unsigned int side = TrkrDefs::getZElement(key);

  // bound check
  if( layer < max_layer && sector < max_phisegment && side < max_zsegment )
  {

    auto ret = m_clusmap[layer][sector][side].insert(std::make_pair(key, newclus));
    if ( !ret.second )
    {
      std::cout << "TrkrClusterContainerv1::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
      exit(1);
    } else {
      ret.first->second->setClusKey( key );
      return ret.first;
    }

  } else {

    std::cout
      << "TrkrClusterContainerv1::addClusterSpecifyKey - out of range access."
      << " layer: " << layer
      << " sector: " << sector
      << " side: " << side
      << std::endl;

    return dummy_map.begin();

  }

}

//_________________________________________________________________
TrkrClusterContainerv1::ConstRange
TrkrClusterContainerv1::getClusters(TrkrDefs::hitsetkey hitsetkey) const
{
  const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  const unsigned int side  = TrkrDefs::getZElement(hitsetkey);
  const unsigned int sector= TrkrDefs::getPhiElement(hitsetkey);

  // bound check
  if( layer < max_layer && sector < max_phisegment && side < max_zsegment )
  {
    return std::make_pair( m_clusmap[layer][sector][side].cbegin(), m_clusmap[layer][sector][side].cend() );

  } else {

    std::cout
      << "TrkrClusterContainerv1::getClusters - out of range access."
      << " layer: " << layer
      << " sector: " << sector
      << " side: " << side
      << std::endl;

    return std::make_pair( dummy_map.cbegin(), dummy_map.cend() );

  }
}

//_________________________________________________________________
TrkrClusterContainerv1::Map*
TrkrClusterContainerv1::getClusterMap(TrkrDefs::hitsetkey hitsetkey)
{
  const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  const unsigned int side  = TrkrDefs::getZElement(hitsetkey);
  const unsigned int sector= TrkrDefs::getPhiElement(hitsetkey);

  // bound check
  if( layer < max_layer && sector < max_phisegment && side < max_zsegment )
  {
    return &m_clusmap[layer][sector][side];
  } else {
    std::cout
      << "TrkrClusterContainerv1::getClusterMap - out of range access."
      << " layer: " << layer
      << " sector: " << sector
      << " side: " << side
      << std::endl;

    return nullptr;
  }
}
  
//_________________________________________________________________
TrkrClusterContainerv1::Iterator
TrkrClusterContainerv1::findOrAddCluster(TrkrDefs::cluskey key)
{
  const unsigned int layer = TrkrDefs::getLayer(key);
  const unsigned int side  = TrkrDefs::getZElement(key);
  const unsigned int sector= TrkrDefs::getPhiElement(key);

  // bound check
  if( layer < max_layer && sector < max_phisegment && side < max_zsegment )
  {

    auto it = m_clusmap[layer][sector][side].lower_bound(key);
    if( it == m_clusmap[layer][sector][side].end() || it->first < key )
    {

      // add new cluster and set its key
      it = m_clusmap[layer][sector][side].insert(it, std::make_pair(key, new TrkrClusterv1()));
      it->second->setClusKey( key );

    }

    return it;
  } else {
    std::cout
      << "TrkrClusterContainerv1::findOrAddCluster - out of range access."
      << " layer: " << layer
      << " sector: " << sector
      << " side: " << side
      << std::endl;

    return dummy_map.begin();
  }

}

//_________________________________________________________________
TrkrCluster* TrkrClusterContainerv1::findCluster(TrkrDefs::cluskey key) const
{
  const unsigned int layer = TrkrDefs::getLayer(key);
  const unsigned int side  = TrkrDefs::getZElement(key);
  const unsigned int sector= TrkrDefs::getPhiElement(key);
  
  // bound check
  if( layer < max_layer && sector < max_phisegment && side < max_zsegment )
  {

    auto it = m_clusmap[layer][sector][side].find(key);
    if (it != m_clusmap[layer][sector][side].end())
    {
      return it->second;
    } else {
      return nullptr;
    }
  
  } else {
    std::cout
      << "TrkrClusterContainerv1::findOrAddCluster - out of range access."
      << " layer: " << layer
      << " sector: " << sector
      << " side: " << side
      << std::endl;
    return nullptr;
  }
}

//_________________________________________________________________
unsigned int TrkrClusterContainerv1::size(void) const
{
  unsigned int size = 0;
  for(unsigned layer = 0;layer < max_layer; layer++)
    for(unsigned phi_segment = 0;phi_segment < max_phisegment;phi_segment++)
    for(unsigned z_segment = 0; z_segment < max_zsegment; z_segment++)
  {
    size += m_clusmap[layer][phi_segment][z_segment].size();
  }
  
  return size;
}
