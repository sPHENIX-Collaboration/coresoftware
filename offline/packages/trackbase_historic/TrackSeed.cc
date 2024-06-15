#include "TrackSeed.h"

TrackSeed::ClusterKeySet dummyKeySet;
TrackSeed::HitIdMap dummyHitIdMap;

TrackSeed::ConstClusterKeyIter TrackSeed::find_cluster_key(TrkrDefs::cluskey /*unused*/) const
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

TrackSeed::ClusterKeyIter TrackSeed::find_cluster_keys(unsigned int /*unused*/)
{
  return dummyKeySet.end();
}

TrackSeed::ClusterKeyIter TrackSeed::end_cluster_keys()
{
  return dummyKeySet.end();
}

TrackSeed::HitIdConstIter TrackSeed::begin_g4hit_id() const
{
  return dummyHitIdMap.end();
}

TrackSeed::HitIdConstIter TrackSeed::find_g4hit_id(int /*unused*/) const
{
  return dummyHitIdMap.end();
}

TrackSeed::HitIdConstIter TrackSeed::end_g4hit_id() const
{
  return dummyHitIdMap.end();
}

TrackSeed::HitIdIter TrackSeed::begin_g4hit_id()
{
  return dummyHitIdMap.end();
}

TrackSeed::HitIdIter TrackSeed::find_g4hit_id(int /*unused*/)
{
  return dummyHitIdMap.end();
}

TrackSeed::HitIdIter TrackSeed::end_g4hit_id()
{
  return dummyHitIdMap.end();
}

const TrackSeed::HitIdMap &TrackSeed::g4hit_ids() const
{
  return dummyHitIdMap;
}
