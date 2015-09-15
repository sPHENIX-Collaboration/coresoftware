#ifndef __SVTXCLUSTERMAP_H__
#define __SVTXCLUSTERMAP_H__

#include "SvtxCluster.h"

#include <phool/PHObject.h>
#include <map>
#include <iostream>

class SvtxClusterMap : public PHObject {
  
public:

  typedef std::map<unsigned int, SvtxCluster*>::const_iterator ConstIter;
  typedef std::map<unsigned int, SvtxCluster*>::iterator            Iter;
  
  SvtxClusterMap();
  SvtxClusterMap(const SvtxClusterMap& clustermap);
  SvtxClusterMap& operator=(const SvtxClusterMap& clustermap);
  virtual ~SvtxClusterMap();

  void identify(std::ostream& os = std::cout) const;
  void Reset();
  int  IsValid() const {return 1;}
  
  bool   empty()                   const {return _map.empty();}
  size_t  size()                   const {return _map.size();}
  size_t count(unsigned int idkey) const {return _map.count(idkey);}
  void   clear()                         {return Reset();}
  
  const SvtxCluster* get(unsigned int idkey) const;
        SvtxCluster* get(unsigned int idkey); 
        SvtxCluster* insert(const SvtxCluster* cluster);
        size_t       erase(unsigned int idkey) {
	  delete _map[idkey]; return _map.erase(idkey);
	}

  ConstIter begin()                   const {return _map.begin();}
  ConstIter  find(unsigned int idkey) const {return _map.find(idkey);}
  ConstIter   end()                   const {return _map.end();}

  Iter begin()                   {return _map.begin();}
  Iter  find(unsigned int idkey) {return _map.find(idkey);}
  Iter   end()                   {return _map.end();}
  
private:
  std::map<unsigned int, SvtxCluster*> _map;
    
  ClassDef(SvtxClusterMap, 1);
};

#endif
