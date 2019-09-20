#ifndef TRACKBASEHISTORIC_SVTXCLUSTERMAPV1_H
#define TRACKBASEHISTORIC_SVTXCLUSTERMAPV1_H

#include "SvtxCluster.h"
#include "SvtxClusterMap.h"

#include <cstddef>
#include <iostream>

class SvtxClusterMap_v1 : public SvtxClusterMap
{
 public:
  SvtxClusterMap_v1();
  SvtxClusterMap_v1(const SvtxClusterMap_v1& clustermap);
  SvtxClusterMap_v1& operator=(const SvtxClusterMap_v1& clustermap);
  virtual ~SvtxClusterMap_v1();

  void identify(std::ostream& os = std::cout) const;
  void Reset();
  int isValid() const { return 1; }
  PHObject* CloneMe() const { return new SvtxClusterMap_v1(*this); }

  bool empty() const { return _map.empty(); }
  size_t size() const { return _map.size(); }
  size_t count(unsigned int idkey) const { return _map.count(idkey); }
  void clear() { return Reset(); }

  const SvtxCluster* get(unsigned int idkey) const;
  SvtxCluster* get(unsigned int idkey);
  SvtxCluster* insert(const SvtxCluster* cluster);
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
  ClusterMap _map;

  ClassDef(SvtxClusterMap_v1, 1);
};

#endif
