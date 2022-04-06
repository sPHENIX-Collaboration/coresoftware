#include "SvtxTrackSeed.h"

SvtxTrackSeed::ClusterKeySet dummyKeySet;


SvtxTrackSeed::ConstClusterKeyIter SvtxTrackSeed::find_cluster_key(TrkrDefs::cluskey /*clusterid*/) const
{
  return dummyKeySet.end();
}

SvtxTrackSeed::ConstClusterKeyIter SvtxTrackSeed::begin_cluster_keys() const
{
  return dummyKeySet.end();
}

SvtxTrackSeed::ConstClusterKeyIter SvtxTrackSeed::end_cluster_keys() const
{
  return dummyKeySet.end();
}


SvtxTrackSeed::ClusterKeyIter SvtxTrackSeed::begin_cluster_keys()
{
  return dummyKeySet.end();
}

SvtxTrackSeed::ClusterKeyIter SvtxTrackSeed::find_cluster_keys(unsigned int /*clusterid*/)
{
  return dummyKeySet.end();
}

SvtxTrackSeed::ClusterKeyIter SvtxTrackSeed::end_cluster_keys()
{
  return dummyKeySet.end();
}
