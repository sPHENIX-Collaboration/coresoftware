#ifndef __SVTXTRACKMAP_H__
#define __SVTXTRACKMAP_H__

#include "SvtxTrack.h"

#include <phool/PHObject.h>
#include <map>

class SvtxTrackMap : public PHObject {
  
public:

  typedef std::map<unsigned int, SvtxTrack*>::const_iterator ConstIter;
  typedef std::map<unsigned int, SvtxTrack*>::iterator            Iter;
  
  SvtxTrackMap();
  SvtxTrackMap(const SvtxTrackMap& trackmap);
  SvtxTrackMap& operator=(const SvtxTrackMap& trackmap);
  virtual ~SvtxTrackMap();

  void identify(std::ostream &os = std::cout) const;
  void Reset();
  int  IsValid() const {return 1;}
  
  bool   empty()                   const {return _map.empty();}
  size_t  size()                   const {return _map.size();}
  size_t count(unsigned int idkey) const {return _map.count(idkey);}
  void   clear()                         {Reset();}
  
  const SvtxTrack* get(unsigned int idkey) const;
        SvtxTrack* get(unsigned int idkey); 
        SvtxTrack* insert(const SvtxTrack* track);
        size_t     erase(unsigned int idkey) {
	  delete _map[idkey]; return _map.erase(idkey);
	}

  ConstIter begin()                   const {return _map.begin();}
  ConstIter  find(unsigned int idkey) const {return _map.find(idkey);}
  ConstIter   end()                   const {return _map.end();}

  Iter begin()                   {return _map.begin();}
  Iter  find(unsigned int idkey) {return _map.find(idkey);}
  Iter   end()                   {return _map.end();}
  
private:
  std::map<unsigned int, SvtxTrack*> _map;
    
  ClassDef(SvtxTrackMap, 1);
};

#endif
