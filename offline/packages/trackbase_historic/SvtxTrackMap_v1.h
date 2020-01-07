#ifndef TRACKBASEHISTORIC_SVTXTRACKMAPV1_H
#define TRACKBASEHISTORIC_SVTXTRACKMAPV1_H

#include "SvtxTrack.h"
#include "SvtxTrackMap.h"

#include <cstddef>        // for size_t
#include <iostream>        // for cout, ostream

class PHObject;

class SvtxTrackMap_v1 : public SvtxTrackMap
{
 public:
  SvtxTrackMap_v1();
  SvtxTrackMap_v1(const SvtxTrackMap_v1& trackmap);
  SvtxTrackMap_v1& operator=(const SvtxTrackMap_v1& trackmap);
  virtual ~SvtxTrackMap_v1();

  void identify(std::ostream& os = std::cout) const;
  void Reset();
  int isValid() const { return 1; }
  PHObject* CloneMe() const { return new SvtxTrackMap_v1(*this); }

  bool empty() const { return _map.empty(); }
  size_t size() const { return _map.size(); }
  size_t count(unsigned int idkey) const { return _map.count(idkey); }
  void clear() { Reset(); }

  const SvtxTrack* get(unsigned int idkey) const;
  SvtxTrack* get(unsigned int idkey);
  SvtxTrack* insert(const SvtxTrack* track);
  size_t erase(unsigned int idkey)
  {
    delete _map[idkey];
    return _map.erase(idkey);
  }

  ConstIter begin() const { return _map.begin(); }
  ConstIter find(unsigned int idkey) const { return _map.find(idkey); }
  ConstIter end() const { return _map.end(); }

  Iter begin() { return _map.begin(); }
  Iter find(unsigned int idkey) { return _map.find(idkey); }
  Iter end() { return _map.end(); }

 private:
  TrackMap _map;

  ClassDef(SvtxTrackMap_v1, 1);
};

#endif
