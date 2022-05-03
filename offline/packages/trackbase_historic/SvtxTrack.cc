#include "SvtxTrack.h"

SvtxTrack::StateMap DummyStateMap;
SvtxTrack::ClusterSet DummyClusterSet;
SvtxTrack::ClusterKeySet DummyClusterKeySet;
SvtxTrack::HitIdMap DummyHitIdMap;

SvtxTrack::HitIdConstIter SvtxTrack::begin_g4hit_id() const
{
  return DummyHitIdMap.end();
}

SvtxTrack::HitIdConstIter SvtxTrack::find_g4hit_id(int /*volume*/) const
{
  return DummyHitIdMap.end();
}

SvtxTrack::HitIdConstIter SvtxTrack::end_g4hit_id() const
{
  return DummyHitIdMap.end();
}


SvtxTrack::HitIdIter SvtxTrack::begin_g4hit_id()
{
  return DummyHitIdMap.end();
}

SvtxTrack::HitIdIter SvtxTrack::find_g4hit_id(int /*volume*/)
{
  return DummyHitIdMap.end();
}

SvtxTrack::HitIdIter SvtxTrack::end_g4hit_id()
{
  return DummyHitIdMap.end();
}


SvtxTrack::ConstStateIter SvtxTrack::begin_states() const
{
  return DummyStateMap.end();
}

SvtxTrack::ConstStateIter SvtxTrack::find_state(float /*pathlength*/)  const
{
  return DummyStateMap.end();
}

SvtxTrack::ConstStateIter SvtxTrack::end_states() const
{
  return DummyStateMap.end();
}

SvtxTrack::StateIter SvtxTrack::begin_states()
{
return DummyStateMap.end();
}

SvtxTrack::StateIter SvtxTrack::find_state(float /*pathlength*/)
{
return DummyStateMap.end();
}

SvtxTrack::StateIter SvtxTrack::end_states()
{
return DummyStateMap.end();
}

SvtxTrack::ConstClusterIter SvtxTrack::begin_clusters() const
{
  return DummyClusterSet.end();
}

SvtxTrack::ConstClusterIter SvtxTrack::find_cluster(unsigned int /*clusterid*/) const
{
  return DummyClusterSet.end();
}

SvtxTrack::ConstClusterIter SvtxTrack::end_clusters() const
{
  return DummyClusterSet.end();
}

SvtxTrack::ClusterIter SvtxTrack::begin_clusters()
{
  return DummyClusterSet.end();
}

SvtxTrack::ClusterIter SvtxTrack::find_cluster(unsigned int /*clusterid*/)
{
  return DummyClusterSet.end();
}

SvtxTrack::ClusterIter SvtxTrack::end_clusters()
{
  return DummyClusterSet.end();
}

SvtxTrack::ConstClusterKeyIter SvtxTrack::find_cluster_key(TrkrDefs::cluskey /*clusterid*/) const
{
  return DummyClusterKeySet.end();
}

SvtxTrack::ConstClusterKeyIter SvtxTrack::begin_cluster_keys() const
{
  return DummyClusterKeySet.end();
}

SvtxTrack::ConstClusterKeyIter SvtxTrack::end_cluster_keys() const
{
  return DummyClusterKeySet.end();
}


SvtxTrack::ClusterKeyIter SvtxTrack::begin_cluster_keys()
{
  return DummyClusterKeySet.end();
}

SvtxTrack::ClusterKeyIter SvtxTrack::find_cluster_keys(unsigned int /*clusterid*/)
{
  return DummyClusterKeySet.end();
}

SvtxTrack::ClusterKeyIter SvtxTrack::end_cluster_keys()
{
  return DummyClusterKeySet.end();
}

const SvtxTrack::HitIdMap &SvtxTrack::g4hit_ids() const
{
  return DummyHitIdMap;
}
