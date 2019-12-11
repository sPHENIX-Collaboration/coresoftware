// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4VERTEX_GLOBALVERTEXMAPV1_H
#define G4VERTEX_GLOBALVERTEXMAPV1_H

#include "GlobalVertexMap.h"

#include "GlobalVertex.h"

#include <cstddef>        // for size_t
#include <iostream>
#include <map>

class GlobalVertexMapv1 : public GlobalVertexMap
{
 public:
  GlobalVertexMapv1();
  virtual ~GlobalVertexMapv1();

  void identify(std::ostream& os = std::cout) const;
  void Reset() { clear(); }
  int isValid() const { return 1; }

  bool empty() const { return _map.empty(); }
  size_t size() const { return _map.size(); }
  size_t count(unsigned int idkey) const { return _map.count(idkey); }
  void clear();

  const GlobalVertex* get(unsigned int idkey) const;
  GlobalVertex* get(unsigned int idkey);
  GlobalVertex* insert(GlobalVertex* vertex);
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
  std::map<unsigned int, GlobalVertex*> _map;

  ClassDef(GlobalVertexMapv1, 1);
};

#endif  // G4VERTEX_GLOBALVERTEXMAPv1_H
