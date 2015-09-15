#ifndef __SVTXHITMAP_H__
#define __SVTXHITMAP_H__

#include "SvtxHit.h"

#include <phool/PHObject.h>
#include <map>
#include <iostream>

class SvtxHitMap : public PHObject {
  
public:

  typedef std::map<unsigned int, SvtxHit*>::const_iterator ConstIter;
  typedef std::map<unsigned int, SvtxHit*>::iterator            Iter;
  
  SvtxHitMap();
  SvtxHitMap(const SvtxHitMap& hitmap);
  SvtxHitMap& operator=(const SvtxHitMap& hitmap);
  virtual ~SvtxHitMap();
  
  void identify(std::ostream& os = std::cout) const;
  void Reset();
  int  IsValid() const {return 1;}
  
  bool   empty()                   const {return _map.empty();}
  size_t  size()                   const {return _map.size();}
  size_t count(unsigned int idkey) const {return _map.count(idkey);}
  void   clear()                         {Reset();}
  
  const SvtxHit* get(unsigned int idkey) const;
        SvtxHit* get(unsigned int idkey); 
        SvtxHit* insert(const SvtxHit *hit);
        size_t   erase(unsigned int idkey) {
	  delete _map[idkey]; return _map.erase(idkey);
	}

  ConstIter begin()                   const {return _map.begin();}
  ConstIter  find(unsigned int idkey) const {return _map.find(idkey);}
  ConstIter   end()                   const {return _map.end();}

  Iter begin()                   {return _map.begin();}
  Iter  find(unsigned int idkey) {return _map.find(idkey);}
  Iter   end()                   {return _map.end();}
  
private:
  std::map<unsigned int, SvtxHit*> _map;
    
  ClassDef(SvtxHitMap, 1);
};

#endif
