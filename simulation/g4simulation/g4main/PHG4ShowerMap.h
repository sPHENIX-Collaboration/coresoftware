#ifndef __PHG4SHOWERMAP_H__
#define __PHG4SHOWERMAP_H__

#include "PHG4Shower.h"

#include <phool/PHObject.h>
#include <map>
#include <iostream>

class PHG4ShowerMap : public PHObject {
  
public:

  typedef std::map<unsigned int, PHG4Shower*>::const_iterator ConstIter;
  typedef std::map<unsigned int, PHG4Shower*>::iterator            Iter;
  
  PHG4ShowerMap();
  virtual ~PHG4ShowerMap();

  void identify(std::ostream &os = std::cout) const;
  void Reset() {clear();}
  int  isValid() const {return 1;}
  
  bool   empty()                   const {return _map.empty();}
  size_t size()                    const {return _map.size();}
  size_t count(unsigned int idkey) const {return _map.count(idkey);}
  void   clear();
  
  const PHG4Shower* get(unsigned int idkey) const;
        PHG4Shower* get(unsigned int idkey); 
        PHG4Shower* insert(PHG4Shower* shower);
        size_t      erase(unsigned int idkey) {
    delete _map[idkey];
    return _map.erase(idkey);
  }

  ConstIter begin()                   const {return _map.begin();}
  ConstIter  find(unsigned int idkey) const {return _map.find(idkey);}
  ConstIter   end()                   const {return _map.end();}

  Iter begin()                   {return _map.begin();}
  Iter  find(unsigned int idkey) {return _map.find(idkey);}
  Iter   end()                   {return _map.end();}
  
private:
  std::map<unsigned int, PHG4Shower*> _map;
    
  ClassDef(PHG4ShowerMap, 1);
};

#endif // __PHG4SHOWERMAP_H__
