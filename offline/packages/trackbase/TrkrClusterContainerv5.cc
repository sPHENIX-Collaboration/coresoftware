/**
 * @file trackbase/TrkrClusterContainerv5.cc
 * @author D. McGlinchey, Hugo Pereira Da Costa
 * @date June 2018
 * @brief Implementation of TrkrClusterContainerv5
 */
#include "TrkrClusterContainerv5.h"
#include <TClonesArray.h>
#include "TrkrCluster.h"
#include "TrkrClusterv4.h"
#include "TrkrClusterv5.h"
#include "TrkrDefs.h"

#include <algorithm>
#include <memory>
#include <string>

using TrkrCluster_impl = TrkrClusterv5;

namespace
{
  TrkrClusterContainer::Map dummy_map;
}

//_________________________________________________________________
void TrkrClusterContainerv5::Reset()
{
  // delete all clusters
  for (auto&& [key, clus_vector] : m_clusmap)
  {
    clus_vector->Clear();
    delete clus_vector;
  }
  // clear the maps
  /* using swap ensures that the memory is properly de-allocated */
  {
    std::map<TrkrDefs::hitsetkey, TClonesArray*> empty;
    m_clusmap.swap(empty);
  }
  // also clear temporary map
  {
    Map empty;
    m_tmpmap.swap(empty);
  }
}

//_________________________________________________________________
void TrkrClusterContainerv5::identify(std::ostream& os) const
{
  os << "-----TrkrClusterContainerv5-----" << std::endl;
  os << "Number of clusters: " << size() << std::endl;

  for (const auto& [hitsetkey, clus_vector] : m_clusmap)
  {
    const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    os << "layer: " << layer << " hitsetkey: " << hitsetkey << std::endl;


    Int_t nClusters = clus_vector->GetEntriesFast();
    for (Int_t i = 0; i < nClusters; i++)
    {
      auto cluster = (TrkrCluster*) (*clus_vector)[i];
      if (cluster)
      {
        cluster->identify(os);
      }
    }
  }

  os << "------------------------------" << std::endl;
}

//_________________________________________________________________
void TrkrClusterContainerv5::removeCluster(TrkrDefs::cluskey key)
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
    const Int_t index = TrkrDefs::getClusIndex(key);

    // compare to vector size
    if (index < clus_vector->GetEntriesFast())
    {
      // delete corresponding element and set to null
      // delete clus_vector[index];
      // clus_vector[index] = nullptr;
      auto removed = clus_vector->RemoveAt(index);
      delete removed;
    }
  }
}

//_________________________________________________________________
void TrkrClusterContainerv5::addClusterSpecifyKey(const TrkrDefs::cluskey key, TrkrCluster* newclus)
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

  // find relevant vector or create one if not found
  // auto clus_vector = std::shared_ptr<TClonesArray>(m_clusmap[hitsetkey]);

  TClonesArray* clus_vector = nullptr;
  auto iter = m_clusmap.find(hitsetkey);
  if (iter != m_clusmap.end())
  {
    clus_vector = iter->second;
  } else {
    // std::string cluster_class_name =
    //     "TrkrClusterv" + std::to_string(m_cluster_version);
    // std::string cluster_class_name = "TrkrCluster";
    TrkrCluster_impl tmp;
    // clus_vector = new TClonesArray(cluster_class_name.c_str());
    clus_vector = new TClonesArray(tmp.GetName());

    m_clusmap[hitsetkey] = clus_vector;
    // std::cout << "created new array for hitsetkey " << hitsetkey << std::endl;
  }

  // get cluster index in vector
  const Int_t index = TrkrDefs::getClusIndex(key);
  Int_t entries = clus_vector->GetEntriesFast();

  // std::cout << "array has size: " << clus_vector->GetEntriesFast() << std::endl;

  TrkrCluster* cls = nullptr;

  // compare index to vector size
  if (index < entries)
  {
    /*
     * if index is already contained in vector, check corresponding element
     * and assign newclus if null
     * print error message and exit otherwise
     */

    auto cluster = (TrkrCluster*) (*clus_vector)[index];
    std::cout << "get null" << std::endl;
    if (!cluster)
    {
      cls = (TrkrCluster_impl*) clus_vector->ConstructedAt(index);
      cls->CopyFrom(static_cast<TrkrCluster_impl*>(newclus));
      // new (clus_vector[index]) TrkrCluster_impl(newclus);
      delete newclus;
    }
    else
    {
      std::cout << "TrkrClusterContainerv5::AddClusterSpecifyKey: duplicate key: " << key << " exiting now" << std::endl;
      exit(1);
    }
  }
  else
  {
    // if index exceeds the vector size, resize cluster to the right size with nullptr, and assign
    // clus_vector.resize(index + 1, nullptr);
    // clus_vector[index] = newclus;
    cls = (TrkrCluster_impl*) clus_vector->ConstructedAt(index);
    cls->CopyFrom(static_cast<TrkrCluster_impl*>(newclus));
    // new (clus_vector[index]) TrkrCluster_impl(newclus);
    delete newclus;
  }
  // std::cout << "cluster valid:" << cls->isValid() << "\n";
}

TrkrClusterContainerv5::ConstRange
TrkrClusterContainerv5::getClusters() const
{
  std::cout << "deprecated function in TrkrClusterContainerv5, user getClusters(TrkrDefs:hitsetkey)"
            << std::endl;
  return std::make_pair(dummy_map.begin(), dummy_map.begin());
}

//_________________________________________________________________
TrkrClusterContainerv5::ConstRange
TrkrClusterContainerv5::getClusters(TrkrDefs::hitsetkey hitsetkey)
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
    for (Int_t index = 0; index < clusters->GetEntriesFast(); ++index)
    {
      const auto& cluster = (TrkrCluster*) (*clusters)[index];
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
TrkrCluster* TrkrClusterContainerv5::findCluster(TrkrDefs::cluskey key) const
{
  // get hitsetkey from cluster
  const TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(key);

  const auto map_iter = m_clusmap.find(hitsetkey);
  if (map_iter != m_clusmap.end())
  {
    // local reference to vector
    const auto& clus_vector = map_iter->second;

    // get cluster position in vector
    const Int_t index = TrkrDefs::getClusIndex(key);

    // compare to vector size
    if (index < clus_vector->GetEntriesFast())
    {
      auto cluster = (TrkrCluster*) (*clus_vector)[index];
      return cluster;
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
TrkrClusterContainer::HitSetKeyList TrkrClusterContainerv5::getHitSetKeys() const
{
  HitSetKeyList out;
  out.reserve(m_clusmap.size());
  std::transform(
      m_clusmap.begin(), m_clusmap.end(), std::back_inserter(out),
      [](const std::pair<TrkrDefs::hitsetkey, TClonesArray*>& pair)
      { return pair.first; });
  return out;
}

//_________________________________________________________________
TrkrClusterContainer::HitSetKeyList TrkrClusterContainerv5::getHitSetKeys(const TrkrDefs::TrkrId trackerid) const
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
      [](const std::pair<TrkrDefs::hitsetkey, TClonesArray*>& pair)
      { return pair.first; });
  return out;
}

//_________________________________________________________________
TrkrClusterContainer::HitSetKeyList TrkrClusterContainerv5::getHitSetKeys(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const
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
      [](const std::pair<TrkrDefs::hitsetkey, TClonesArray*>& pair)
      { return pair.first; });
  return out;
}

//_________________________________________________________________
unsigned int TrkrClusterContainerv5::size() const
{
  unsigned int size = 0;
  for (const auto& [hitsetkey, clus_vector] : m_clusmap)
  {
    size += (unsigned int) clus_vector->GetEntries();
  }
  return size;
}
