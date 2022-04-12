#include "TrackSeed.h"

TrackSeed::ClusterKeySet dummyKeySet;


TrackSeed::ConstClusterKeyIter TrackSeed::find_cluster_key(TrkrDefs::cluskey) const
{
  return dummyKeySet.end();
}

TrackSeed::ConstClusterKeyIter TrackSeed::begin_cluster_keys() const
{
  return dummyKeySet.end();
}

TrackSeed::ConstClusterKeyIter TrackSeed::end_cluster_keys() const
{
  return dummyKeySet.end();
}


TrackSeed::ClusterKeyIter TrackSeed::begin_cluster_keys()
{
  return dummyKeySet.end();
}

TrackSeed::ClusterKeyIter TrackSeed::find_cluster_keys(unsigned int)
{
  return dummyKeySet.end();
}

TrackSeed::ClusterKeyIter TrackSeed::end_cluster_keys()
{
  return dummyKeySet.end();
}
