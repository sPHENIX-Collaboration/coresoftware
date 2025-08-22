#include "PhotonClusterContainer.h"

#include "PhotonCluster.h"
#include "PhotonClusterv1.h" // for concrete type id access (set_id via RawCluster side not available here)

#include <iostream>

PhotonClusterContainer::ConstRange PhotonClusterContainer::getClusters() const
{
  return std::make_pair(m_clusters.begin(), m_clusters.end());
}

PhotonClusterContainer::Range PhotonClusterContainer::getClusters()
{
  return std::make_pair(m_clusters.begin(), m_clusters.end());
}

PhotonClusterContainer::ConstIterator PhotonClusterContainer::AddCluster(PhotonCluster* clus)
{
  unsigned int key = m_clusters.size();
  while (m_clusters.find(key) != m_clusters.end())
  {
    key++;
  }
  // We cannot set id directly because PhotonCluster interface has no set_id; assume dynamic_cast to RawClusterv1 if needed externally.
  m_clusters[key] = clus;
  return m_clusters.find(key);
}

PhotonCluster* PhotonClusterContainer::getCluster(const unsigned int key)
{
  ConstIterator it = m_clusters.find(key);
  if (it != m_clusters.end()) return it->second;
  return nullptr;
}

const PhotonCluster* PhotonClusterContainer::getCluster(const unsigned int key) const
{
  ConstIterator it = m_clusters.find(key);
  if (it != m_clusters.end()) return it->second;
  return nullptr;
}

int PhotonClusterContainer::isValid() const { return (!m_clusters.empty()); }

void PhotonClusterContainer::Reset()
{
  while (!m_clusters.empty())
  {
    delete m_clusters.begin()->second;
    m_clusters.erase(m_clusters.begin());
  }
}

void PhotonClusterContainer::identify(std::ostream& os) const
{
  os << "PhotonClusterContainer, number of clusters: " << size() << std::endl;
}
