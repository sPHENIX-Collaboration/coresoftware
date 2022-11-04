/**
 * @file trackbase/TrkrClusterContainerv3.cc
 * @author D. McGlinchey, Hugo Pereira Da Costa
 * @date June 2018
 * @brief Implementation of TrkrClusterContainerv3
 */
#include "TrkrClusterContainerv3.h"
#include "TrkrCluster.h"
#include "TrkrDefs.h"

#include <algorithm>

namespace
{
  TrkrClusterContainer::Map dummy_map;
}

//_________________________________________________________________
void TrkrClusterContainerv3::Reset()
{
  // delete all clusters
  for (auto&& [key, map] : m_clusmap)
  {
    for (auto&& [cluskey, cluster] : map)
    {
      delete cluster;
    }
  }

  // clear the maps
  /* using swap ensures that the memory is properly de-allocated */
  std::map<TrkrDefs::hitsetkey, Map> empty;
  m_clusmap.swap(empty);
}

//_________________________________________________________________
void TrkrClusterContainerv3::identify(std::ostream& os) const
{
  os << "-----TrkrClusterContainerv3-----" << std::endl;
  os << "Number of clusters: " << size() << std::endl;

  for (const auto& [hitsetkey, map] : m_clusmap)
  {
    const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    os << "layer: " << layer << " hitsetkey: " << hitsetkey << std::endl;
    for (const auto& [cluskey, cluster] : map)
    {
      os << "clus key " << cluskey << " layer " << TrkrDefs::getLayer(cluskey) << std::endl;
      cluster->identify();
    }
  }

  os << "------------------------------" << std::endl;
}

//_________________________________________________________________
void TrkrClusterContainerv3::removeCluster(TrkrDefs::cluskey key)
{
  // get hitset key from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

  // find relevant cluster map if any and remove corresponding cluster
  auto iter = m_clusmap.find(hitsetkey);
  if (iter != m_clusmap.end())
  {
    TrkrCluster* clus = findCluster(key);
    delete clus;
    iter->second.erase(key);
  }
}

//_________________________________________________________________
void TrkrClusterContainerv3::addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

  // find relevant cluster map or create one if not found
  Map& map = m_clusmap[hitsetkey];
  const auto [iter, success] = map.insert(std::make_pair(key, newclus));
  if (!success)
  {
    std::cout << "TrkrClusterContainerv3::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
    exit(1);
  }
}

//_________________________________________________________________
TrkrClusterContainerv3::ConstRange
TrkrClusterContainerv3::getClusters(TrkrDefs::hitsetkey hitsetkey)
{
  // find relevant association map
  const auto iter = m_clusmap.find(hitsetkey);
  if (iter != m_clusmap.end())
  {
    return std::make_pair(iter->second.cbegin(), iter->second.cend());
  }
  else
  {
    return std::make_pair(dummy_map.cbegin(), dummy_map.cend());
  }
}

//_________________________________________________________________
TrkrCluster* TrkrClusterContainerv3::findCluster(TrkrDefs::cluskey key) const
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

  const auto map_iter = m_clusmap.find(hitsetkey);
  if (map_iter != m_clusmap.end())
  {
    const auto clus_iter = map_iter->second.find(key);
    if (clus_iter != map_iter->second.end())
    {
      return clus_iter->second;
    }
    else
    {
      return nullptr;
    }
  }
  else
  {
    return nullptr;
  }
}

//_________________________________________________________________
TrkrClusterContainer::HitSetKeyList TrkrClusterContainerv3::getHitSetKeys() const
{
  HitSetKeyList out;
  out.reserve(m_clusmap.size());
  std::transform(
      m_clusmap.begin(), m_clusmap.end(), std::back_inserter(out),
      [](const std::pair<TrkrDefs::hitsetkey, Map>& pair)
      { return pair.first; });
  return out;
}

//_________________________________________________________________
TrkrClusterContainer::HitSetKeyList TrkrClusterContainerv3::getHitSetKeys(const TrkrDefs::TrkrId trackerid) const
{
  /* copy the logic from TrkrHitSetContainerv1::getHitSets */
  const TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid);
  const TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid);

  // get relevant range in map
  const auto begin = m_clusmap.lower_bound(keylo);
  const auto end = m_clusmap.upper_bound(keyhi);

  // transform to a vector
  HitSetKeyList out;
  out.reserve(m_clusmap.size());
  std::transform(
      begin, end, std::back_inserter(out),
      [](const std::pair<TrkrDefs::hitsetkey, Map>& pair)
      { return pair.first; });
  return out;
}

//_________________________________________________________________
TrkrClusterContainer::HitSetKeyList TrkrClusterContainerv3::getHitSetKeys(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const
{
  /* copy the logic from TrkrHitSetContainerv1::getHitSets */
  TrkrDefs::hitsetkey keylo = TrkrDefs::getHitSetKeyLo(trackerid, layer);
  TrkrDefs::hitsetkey keyhi = TrkrDefs::getHitSetKeyHi(trackerid, layer);

  // get relevant range in map
  const auto begin = m_clusmap.lower_bound(keylo);
  const auto end = m_clusmap.upper_bound(keyhi);

  // transform to a vector
  HitSetKeyList out;
  out.reserve(m_clusmap.size());
  std::transform(
      begin, end, std::back_inserter(out),
      [](const std::pair<TrkrDefs::hitsetkey, Map>& pair)
      { return pair.first; });
  return out;
}

//_________________________________________________________________
unsigned int TrkrClusterContainerv3::size() const
{
  unsigned int size = 0;
  for (const auto& map_pair : m_clusmap)
  {
    size += map_pair.second.size();
  }

  return size;
}
