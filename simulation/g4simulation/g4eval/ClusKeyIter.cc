#include "ClusKeyIter.h"

#include <trackbase_historic/SvtxTrack.h>

ClusKeyIter::ClusKeyIter(SvtxTrack* _track)
  : track{_track}
  , in_silicon{_track->get_silicon_seed() != nullptr}
  , has_tpc{_track->get_tpc_seed() != nullptr}
  , no_data{!in_silicon && !has_tpc}
{
}

ClusKeyIter ClusKeyIter::begin() const
{
  ClusKeyIter iter0{track};
  if (iter0.no_data)
  {
    return iter0;
  }
  if (iter0.in_silicon)
  {
    iter0.iter = track->get_silicon_seed()->begin_cluster_keys();
    iter0.iter_end_silicon = track->get_silicon_seed()->end_cluster_keys();
  }
  else if (has_tpc)
  {
    iter0.iter = track->get_tpc_seed()->begin_cluster_keys();
  }
  return iter0;
}

ClusKeyIter ClusKeyIter::end() const
{
  ClusKeyIter iter0{track};
  if (iter0.no_data)
  {
    return iter0;
  }
  if (has_tpc)
  {
    iter0.iter = track->get_tpc_seed()->end_cluster_keys();
  }
  else if (in_silicon)
  {
    iter0.iter = track->get_silicon_seed()->end_cluster_keys();
  }
  return iter0;
}

void ClusKeyIter::operator++()
{
  if (no_data)
  {
    return;
  }
  ++iter;
  if (in_silicon && has_tpc && iter == iter_end_silicon)
  {
    in_silicon = false;
    iter = track->get_tpc_seed()->begin_cluster_keys();
  }
}

bool ClusKeyIter::operator!=(const ClusKeyIter& rhs) const
{
  if (no_data)
  {
    return false;
  }
  return iter != rhs.iter;
}

TrkrDefs::cluskey ClusKeyIter::operator*() const
{
  return *iter;
}
