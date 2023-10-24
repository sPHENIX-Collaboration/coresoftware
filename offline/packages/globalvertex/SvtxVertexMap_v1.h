#ifndef TRACKBASEHISTORIC_SVTXVERTEXMAPV1_H
#define TRACKBASEHISTORIC_SVTXVERTEXMAPV1_H

#include "SvtxVertex.h"
#include "SvtxVertexMap.h"

#include <cstddef>         // for size_t
#include <iostream>
#include <map>

class PHObject;

class SvtxVertexMap_v1 : public SvtxVertexMap
{
 public:
  SvtxVertexMap_v1();
  SvtxVertexMap_v1(const SvtxVertexMap_v1& vertexmap);
  SvtxVertexMap_v1& operator=(const SvtxVertexMap_v1& vertexmap);
  ~SvtxVertexMap_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override;
  int isValid() const override { return 1; }
  PHObject* CloneMe() const override { return new SvtxVertexMap_v1(*this); }

  bool empty() const override { return _map.empty(); }
  size_t size() const override { return _map.size(); }
  size_t count(unsigned int idkey) const override { return _map.count(idkey); }
  void clear() override { Reset(); }

  const SvtxVertex* get(unsigned int idkey) const override;
  SvtxVertex* get(unsigned int idkey) override;

  //! Add vertex to container. Note the container takes ownership
  SvtxVertex* insert(SvtxVertex* vertex) override;
  //! legacy interface. Add vertex to container. Note the container does not take ownership
  SvtxVertex* insert_clone(const SvtxVertex* vertex) override;
  size_t erase(unsigned int idkey) override
  {
    delete _map[idkey];
    return _map.erase(idkey);
  }

  ConstIter begin() const override { return _map.begin(); }
  ConstIter find(unsigned int idkey) const override { return _map.find(idkey); }
  ConstIter end() const override { return _map.end(); }

  Iter begin() override { return _map.begin(); }
  Iter find(unsigned int idkey) override { return _map.find(idkey); }
  Iter end() override { return _map.end(); }

 private:
  std::map<unsigned int, SvtxVertex*> _map;

  ClassDefOverride(SvtxVertexMap_v1, 1);
};

#endif  // __SVTXVERTEXMAP_V1_H__
