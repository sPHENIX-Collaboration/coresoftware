// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_CALOVERTEXMAPV1_H
#define GLOBALVERTEX_CALOVERTEXMAPV1_H

#include "CaloVertexMap.h"

#include "CaloVertex.h"

#include <cstddef>  // for size_t
#include <iostream>
#include <map>

class CaloVertexMapv1 : public CaloVertexMap
{
 public:
  CaloVertexMapv1() = default;
  ~CaloVertexMapv1() override;

  void identify(std::ostream& os = std::cout) const override;
  void Reset() override { clear(); }
  int isValid() const override { return 1; }

  bool empty() const override { return _map.empty(); }
  size_t size() const override { return _map.size(); }
  size_t count(unsigned int idkey) const override { return _map.count(idkey); }
  void clear() override;

  const CaloVertex* get(unsigned int idkey) const override;
  CaloVertex* get(unsigned int idkey) override;
  CaloVertex* insert(CaloVertex* vertex) override;
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
  std::map<unsigned int, CaloVertex*> _map;

  ClassDefOverride(CaloVertexMapv1, 1);
};

#endif  // G4CALO_CALOVERTEXMAPV1_H
