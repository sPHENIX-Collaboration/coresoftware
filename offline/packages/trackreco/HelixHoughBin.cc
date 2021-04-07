#include "HelixHoughBin.h"

HelixHoughBin::ClusterSet HelixHoughBinClusterSet;

HelixHoughBin::ConstClusterIter HelixHoughBin::begin_clusters() const
{
  return HelixHoughBinClusterSet.end();
}

HelixHoughBin::ConstClusterIter HelixHoughBin::find_cluster(unsigned int cluster_id) const
{
  return HelixHoughBinClusterSet.end();
}
HelixHoughBin::ConstClusterIter HelixHoughBin::end_clusters() const
{
  return HelixHoughBinClusterSet.end();
}

HelixHoughBin::ClusterIter HelixHoughBin::begin_clusters()
{
  return HelixHoughBinClusterSet.end();
}

HelixHoughBin::ClusterIter HelixHoughBin::find_cluster(unsigned int cluster_id)
{
  return HelixHoughBinClusterSet.end();
}

HelixHoughBin::ClusterIter HelixHoughBin::end_clusters()
{
  return HelixHoughBinClusterSet.end();
}
