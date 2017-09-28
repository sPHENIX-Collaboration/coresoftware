#ifndef __PHHEPMCGENEVENTMAP_H__
#define __PHHEPMCGENEVENTMAP_H__

#include "PHHepMCGenEvent.h"

#include <phool/PHObject.h>
#include <map>

class PHHepMCGenEventMap : public PHObject {
  
public:

  typedef std::map<unsigned int, PHHepMCGenEvent*> HepMCGenEventMap;
  typedef std::map<unsigned int, PHHepMCGenEvent*>::const_iterator ConstIter;
  typedef std::map<unsigned int, PHHepMCGenEvent*>::iterator            Iter;  
  
  PHHepMCGenEventMap();
  PHHepMCGenEventMap(const PHHepMCGenEventMap& eventmap);
  PHHepMCGenEventMap& operator=(const PHHepMCGenEventMap& eventmap);
  virtual ~PHHepMCGenEventMap();

  void identify(std::ostream &os = std::cout) const;
  void Reset();
  int  isValid() const {return 1;}
  PHHepMCGenEventMap* Clone() const {return new PHHepMCGenEventMap(*this);}
  
  bool   empty()                   const {return _map.empty();}
  size_t  size()                   const {return _map.size();}
  size_t count(unsigned int idkey) const {return _map.count(idkey);}
  void   clear()                         {Reset();}
  
  const PHHepMCGenEvent* get(unsigned int idkey) const;
        PHHepMCGenEvent* get(unsigned int idkey); 
        PHHepMCGenEvent* insert(const PHHepMCGenEvent* event);
        size_t      erase(unsigned int idkey) {
	  delete _map[idkey]; return _map.erase(idkey);
	}

  ConstIter begin()                   const {return _map.begin();}
  ConstIter  find(unsigned int idkey) const {return _map.find(idkey);}
  ConstIter   end()                   const {return _map.end();}

  Iter begin()                   {return _map.begin();}
  Iter  find(unsigned int idkey) {return _map.find(idkey);}
  Iter   end()                   {return _map.end();}
  
private:
  HepMCGenEventMap _map;
    
  ClassDef(PHHepMCGenEventMap, 1);
};

#endif
