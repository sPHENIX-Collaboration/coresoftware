#ifndef __SVTXVERTEXMAP_H__
#define __SVTXVERTEXMAP_H__

#include "SvtxVertex.h"

#include <phool/PHObject.h>
#include <map>
#include <iostream>

class SvtxVertexMap : public PHObject {
  
public:

  typedef std::map<unsigned int, SvtxVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, SvtxVertex*>::iterator            Iter;
  
  SvtxVertexMap();
  virtual ~SvtxVertexMap();

  void identify(std::ostream &os = std::cout) const;
  void Reset() {clear();}
  int  IsValid() const {return 1;}
  
  bool   empty()                   const {return _map.empty();}
  size_t  size()                   const {return _map.size();}
  size_t count(unsigned int idkey) const {return _map.count(idkey);}
  void   clear()                         {return _map.clear();}
  
  const SvtxVertex* get(unsigned int idkey) const;
        SvtxVertex* get(unsigned int idkey); 
        SvtxVertex* insert(const SvtxVertex* vertex);
        size_t      erase(unsigned int idkey) {return _map.erase(idkey);}

  ConstIter begin()                   const {return _map.begin();}
  ConstIter  find(unsigned int idkey) const {return _map.find(idkey);}
  ConstIter   end()                   const {return _map.end();}

  Iter begin()                   {return _map.begin();}
  Iter  find(unsigned int idkey) {return _map.find(idkey);}
  Iter   end()                   {return _map.end();}
  
private:
  std::map<unsigned int, SvtxVertex*> _map;
    
  ClassDef(SvtxVertexMap, 1);
};

#endif // __SVTXVERTEXLIST_H__
