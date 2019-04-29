#ifndef __SVTXVERTEXMAP_V1_H__
#define __SVTXVERTEXMAP_V1_H__

#include "SvtxVertex.h"
#include "SvtxVertexMap.h"

#include <phool/PHObject.h>
#include <iostream>
#include <map>

class SvtxVertexMap_v1 : public SvtxVertexMap
{
 public:
  SvtxVertexMap_v1();
  SvtxVertexMap_v1(const SvtxVertexMap_v1& vertexmap);
  SvtxVertexMap_v1& operator=(const SvtxVertexMap_v1& vertexmap);
  virtual ~SvtxVertexMap_v1();

  void identify(std::ostream& os = std::cout) const;
  void Reset();
  int isValid() const { return 1; }
  SvtxVertexMap* Clone() const { return new SvtxVertexMap_v1(*this); }

  bool empty() const { return _map.empty(); }
  size_t size() const { return _map.size(); }
  size_t count(unsigned int idkey) const { return _map.count(idkey); }
  void clear() { Reset(); }

  const SvtxVertex* get(unsigned int idkey) const;
  SvtxVertex* get(unsigned int idkey);

  //! Add vertex to container. Note the container takes ownership
  SvtxVertex* insert(SvtxVertex* vertex);
  //! legacy interface. Add vertex to container. Note the container does not take ownership
  SvtxVertex* insert_clone(const SvtxVertex* vertex);
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
  std::map<unsigned int, SvtxVertex*> _map;

  ClassDef(SvtxVertexMap_v1, 1);
};

#endif  // __SVTXVERTEXMAP_V1_H__
