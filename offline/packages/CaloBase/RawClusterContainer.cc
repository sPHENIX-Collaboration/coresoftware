#include "RawClusterContainer.h"

#include "RawCluster.h"

#include <cstdlib>
#include <iostream>

using namespace std;

RawClusterContainer::ConstRange
RawClusterContainer::getClusters() const
{
  return make_pair(_clusters.begin(), _clusters.end());
}

RawClusterContainer::Range
RawClusterContainer::getClusters()
{
  return make_pair(_clusters.begin(), _clusters.end());
}

RawClusterContainer::ConstIterator
RawClusterContainer::AddCluster(RawCluster* rawcluster)
{
  unsigned int key = _clusters.size();
  // just to be safe in case someone deleted a cluster and key is
  // a valid index of a cluster, increment key until there is no cluster
  // in our map
  while (_clusters.find(key) != _clusters.end())
  {
    key++;
  }
  rawcluster->set_id(key);
  _clusters[key] = rawcluster;
  return _clusters.find(key);
}

RawCluster*
RawClusterContainer::getCluster(const unsigned int key)
{
  ConstIterator it = _clusters.find(key);
  if (it != _clusters.end())
  {
    return it->second;
  }
  return nullptr;
}

const RawCluster*
RawClusterContainer::getCluster(const unsigned int key) const
{
  ConstIterator it = _clusters.find(key);
  if (it != _clusters.end())
  {
    return it->second;
  }
  return nullptr;
}

int RawClusterContainer::isValid() const
{
  return (!_clusters.empty());
}

void RawClusterContainer::Reset()
{
  while (_clusters.begin() != _clusters.end())
  {
    delete _clusters.begin()->second;
    _clusters.erase(_clusters.begin());
  }
}

void RawClusterContainer::identify(std::ostream& os) const
{
  os << "RawClusterContainer, number of clusters: " << size() << std::endl;
}

double
RawClusterContainer::getTotalEdep() const
{
  double totalenergy = 0;
  ConstIterator iter;
  for (iter = _clusters.begin(); iter != _clusters.end(); ++iter)
  {
    totalenergy += iter->second->get_energy();
  }
  return totalenergy;
}
