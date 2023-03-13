/**
 * @file trackbase/TrkrClusterContainerv4.cc
 * @author D. McGlinchey, Hugo Pereira Da Costa
 * @date June 2018
 * @brief Implementation of TrkrClusterContainerv4
 */
#include "TrkrClusterContainerv4.h"
#include "TrkrCluster.h"
#include "TrkrDefs.h"

#include <algorithm>

namespace
{
  TrkrClusterContainer::Map dummy_map;
}

//_________________________________________________________________
void TrkrClusterContainerv4::Reset()
{
  // delete all clusters
  for (auto&& [key, clus_vector] : m_clusmap)
  {
    for (auto&& cluster : clus_vector)
    {
      delete cluster;
    }
  }

  // clear the maps
  /* using swap ensures that the memory is properly de-allocated */
  {
    std::map<TrkrDefs::hitsetkey, Vector> empty;
    m_clusmap.swap(empty);
  }

  // also clear temporary map
  {
    Map empty;
    m_tmpmap.swap(empty);
  }
}

//_________________________________________________________________
void TrkrClusterContainerv4::identify(std::ostream& os) const
{
  os << "-----TrkrClusterContainerv4-----" << std::endl;
  os << "Number of clusters: " << size() << std::endl;

  for (const auto& [hitsetkey, clus_vector] : m_clusmap)
  {
    const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    os << "layer: " << layer << " hitsetkey: " << hitsetkey << std::endl;

    for (const auto& cluster : clus_vector)
    {
      if (cluster)
      {
        cluster->identify(os);
      }
    }
  }

  os << "------------------------------" << std::endl;
}

//_________________________________________________________________
void TrkrClusterContainerv4::removeCluster(TrkrDefs::cluskey key)
{
  // get hitset key from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

  // find relevant cluster map if any and remove corresponding cluster
  auto iter = m_clusmap.find(hitsetkey);
  if (iter != m_clusmap.end())
  {
    // local reference to the vector
    auto& clus_vector = iter->second;

    // cluster index in vector
    const auto index = TrkrDefs::getClusIndex(key);

    // compare to vector size
    if (index < clus_vector.size())
    {
      // delete corresponding element and set to null
      delete clus_vector[index];
      clus_vector[index] = nullptr;
    }
  }
}

//_________________________________________________________________
void TrkrClusterContainerv4::addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

  // find relevant vector or create one if not found
  auto& clus_vector = m_clusmap[hitsetkey];

  // get cluster index in vector
  const auto index = TrkrDefs::getClusIndex(key);

  // compare index to vector size
  if (index < clus_vector.size())
  {
    /*
     * if index is already contained in vector, check corresponding element
     * and assign newclus if null
     * print error message and exit otherwise
     */
    if (!clus_vector[index])
    {
      clus_vector[index] = newclus;
    }
    else
    {
      std::cout << "TrkrClusterContainerv4::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
      exit(1);
    }
  }
  else if (index == clus_vector.size())
  {
    // if index matches the vector size, just push back the new cluster
    clus_vector.push_back(newclus);
  }
  else
  {
    // if index exceeds the vector size, resize cluster to the right size with nullptr, and assign
    clus_vector.resize(index + 1, nullptr);
    clus_vector[index] = newclus;
  }
}

TrkrClusterContainerv4::ConstRange
TrkrClusterContainerv4::getClusters() const
{
  std::cout << "deprecated function in TrkrClusterContainerv4, user getClusters(TrkrDefs:hitsetkey)"
            << std::endl;
  return std::make_pair(dummy_map.begin(), dummy_map.begin());
}

//_________________________________________________________________
TrkrClusterContainerv4::ConstRange
TrkrClusterContainerv4::getClusters(TrkrDefs::hitsetkey hitsetkey)
{
  // clear temporary map
  {
    Map empty;
    m_tmpmap.swap(empty);
  }

  // find relevant vector
  const auto iter = m_clusmap.find(hitsetkey);
  if (iter != m_clusmap.end())
  {
    // copy content in temporary map
    const auto& clusters = iter->second;
    for (size_t index = 0; index < clusters.size(); ++index)
    {
      const auto& cluster = clusters[index];
      if (cluster)
      {
        // generate cluster key from hitset and index
        const auto ckey = TrkrDefs::genClusKey(hitsetkey, index);

        // insert in map
        m_tmpmap.insert(m_tmpmap.end(), std::make_pair(ckey, cluster));
      }
    }
  }

  // return temporary map range
  return std::make_pair(m_tmpmap.cbegin(), m_tmpmap.cend());
}

//_________________________________________________________________
TrkrCluster* TrkrClusterContainerv4::findCluster(TrkrDefs::cluskey key) const
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

  const auto map_iter = m_clusmap.find(hitsetkey);
  if (map_iter != m_clusmap.end())
  {
    // local reference to vector
    const auto& clus_vector = map_iter->second;

    // get cluster position in vector
    const auto index = TrkrDefs::getClusIndex(key);

    // compare to vector size
    if (index < clus_vector.size())
    {
      return clus_vector[index];
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
TrkrClusterContainer::HitSetKeyList TrkrClusterContainerv4::getHitSetKeys() const
{
  HitSetKeyList out;
  out.reserve(m_clusmap.size());
  std::transform(
      m_clusmap.begin(), m_clusmap.end(), std::back_inserter(out),
      [](const std::pair<TrkrDefs::hitsetkey, Vector>& pair)
      { return pair.first; });
  return out;
}

//_________________________________________________________________
TrkrClusterContainer::HitSetKeyList TrkrClusterContainerv4::getHitSetKeys(const TrkrDefs::TrkrId trackerid) const
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
      [](const std::pair<TrkrDefs::hitsetkey, Vector>& pair)
      { return pair.first; });
  return out;
}

//_________________________________________________________________
TrkrClusterContainer::HitSetKeyList TrkrClusterContainerv4::getHitSetKeys(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const
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
      [](const std::pair<TrkrDefs::hitsetkey, Vector>& pair)
      { return pair.first; });
  return out;
}

//_________________________________________________________________
unsigned int TrkrClusterContainerv4::size() const
{
  unsigned int size = 0;
  for (const auto& [hitsetkey, clus_vector] : m_clusmap)
  {
    size += std::count_if(clus_vector.begin(), clus_vector.end(), [](TrkrCluster* cluster)
                          { return cluster; });
  }
  return size;
}
