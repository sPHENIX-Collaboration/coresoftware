#ifndef GLOBALVERTEX_TRUTHVERTEXMAPV1_H
#define GLOBALVERTEX_TRUTHVERTEXMAPV1_H

#include "TruthVertex.h"
#include "TruthVertexMap.h"

#include <iostream>
#include <map>

class TruthVertexMap_v1 : public TruthVertexMap
{
 public:
  TruthVertexMap_v1() = default;
  ~TruthVertexMap_v1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { clear(); }
  int isValid() const override { return !_map.empty(); }

  bool empty() const override { return _map.empty(); }
  size_t size() const override { return _map.size(); }
  size_t count(unsigned int idkey) const override { return _map.count(idkey); }
  void clear() override;

  const TruthVertex* get(unsigned int idkey) const override;
  TruthVertex* get(unsigned int idkey) override;
  TruthVertex* insert(TruthVertex* vertex) override;
  size_t erase(unsigned int idkey) override;

  ConstIter begin() const override { return _map.begin(); }
  ConstIter find(unsigned int idkey) const override { return _map.find(idkey); }
  ConstIter end() const override { return _map.end(); }

  Iter begin() override { return _map.begin(); }
  Iter find(unsigned int idkey) override { return _map.find(idkey); }
  Iter end() override { return _map.end(); }

 private:
  std::map<unsigned int, TruthVertex*> _map;

  ClassDefOverride(TruthVertexMap_v1, 1);
};

#endif
