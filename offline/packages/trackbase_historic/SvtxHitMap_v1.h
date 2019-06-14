#ifndef TRACKBASEHISTORIC_SVTXHITMAPV1_H
#define TRACKBASEHISTORIC_SVTXHITMAPV1_H

#include "SvtxHit.h"
#include "SvtxHitMap.h"

#include <cstddef>
#include <iostream>

class SvtxHitMap_v1 : public SvtxHitMap
{
 public:
  SvtxHitMap_v1();
  SvtxHitMap_v1(const SvtxHitMap_v1& hitmap);
  SvtxHitMap_v1& operator=(const SvtxHitMap_v1& hitmap);
  virtual ~SvtxHitMap_v1();

  void identify(std::ostream& os = std::cout) const;
  void Reset();
  int isValid() const { return 1; }
  SvtxHitMap* clone() const { return new SvtxHitMap_v1(*this); }

  bool empty() const { return _map.empty(); }
  size_t size() const { return _map.size(); }
  size_t count(unsigned int idkey) const { return _map.count(idkey); }
  void clear() { Reset(); }

  const SvtxHit* get(unsigned int idkey) const;
  SvtxHit* get(unsigned int idkey);
  SvtxHit* insert(const SvtxHit* hit);
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
  HitMap _map;

  ClassDef(SvtxHitMap_v1, 1);
};

#endif
