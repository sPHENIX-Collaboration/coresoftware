#ifndef CLUSKEYITER__H
#define CLUSKEYITER__H

// an iterator to loop over all the TrkrClusters for a given track
#include <set>
#include <trackbase/TrkrDefs.h>

class SvtxTrack;

struct ClusKeyIter {
  typedef std::set<TrkrDefs::cluskey> ClusterKeySet;
  typedef ClusterKeySet::iterator ClusterKeyIter;

  ClusKeyIter(SvtxTrack* _track);
  // data
  SvtxTrack* track;
  bool in_silicon;
  bool has_tpc;
  bool no_data; // neither a tpc nor a silicon seed
  ClusterKeyIter iter             { };
  ClusterKeyIter iter_end_silicon { };

  ClusKeyIter begin();
  ClusKeyIter end();

  void operator++();
  TrkrDefs::cluskey operator*();
  bool operator!=(const ClusKeyIter& rhs);
};



#endif
