#ifndef CLUSKEYITER__H
#define CLUSKEYITER__H

// an iterator to loop over all the TrkrClusters for a given track
#include <trackbase/TrkrDefs.h>
#include <set>

class SvtxTrack;

struct ClusKeyIter
{
  typedef std::set<TrkrDefs::cluskey> ClusterKeySet;
  typedef ClusterKeySet::iterator ClusterKeyIter;

  ClusKeyIter(SvtxTrack* _track);
  // data
  SvtxTrack* track;
  bool in_silicon;
  bool has_tpc;
  bool no_data;  // neither a tpc nor a silicon seed
  ClusterKeyIter iter{};
  ClusterKeyIter iter_end_silicon{};

  ClusKeyIter begin() const;
  ClusKeyIter end() const;

  void operator++();
  TrkrDefs::cluskey operator*() const;
  bool operator!=(const ClusKeyIter& rhs) const;
};

#endif
