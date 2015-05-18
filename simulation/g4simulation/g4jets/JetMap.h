#ifndef __JETMAP_H__
#define __JETMAP_H__

#include "Jet.h"

#include <phool/PHObject.h>
#include <set>
#include <map>
#include <cmath>

class JetMap : public PHObject {
  
public:

  // jet object iterators
  typedef std::map<unsigned int, Jet*>  typ_JetMap;
  typedef typ_JetMap::const_iterator    ConstIter;
  typedef typ_JetMap::iterator          Iter;

  // source identifier iterators
  typedef std::set<Jet::SRC>::const_iterator ConstSrcIter;
  typedef std::set<Jet::SRC>::iterator            SrcIter;
  
  JetMap();
  virtual ~JetMap();

  void identify(std::ostream &os = std::cout) const;
  void Reset();
  int  isValid() const {return 1;}

  // map content info
  
  void set_algo(Jet::ALGO algo) {_algo = algo;}
  Jet::ALGO get_algo() const    {return _algo;}
  
  void  set_par(float par) {_par = par;}
  float get_par() const    {return _par;}

  // set access to list of source identifiers
  
  bool empty_src() const        {return _src.empty();}
  void insert_src(Jet::SRC src) {_src.insert(src);}
 
  ConstSrcIter begin_src()             const {return _src.begin();}
  ConstSrcIter  find_src(Jet::SRC src) const {return _src.find(src);}
  ConstSrcIter   end_src()             const {return _src.end();}

  SrcIter begin_src()             {return _src.begin();}
  SrcIter  find_src(Jet::SRC src) {return _src.find(src);}
  SrcIter   end_src()             {return _src.end();}
  
  // map access to jets
 
  bool   empty()                   const {return _map.empty();}
  size_t  size()                   const {return _map.size();}
  size_t count(unsigned int idkey) const {return _map.count(idkey);}
  void   clear()                         {Reset();}
  
  const Jet*   get(unsigned int idkey) const;
        Jet*   get(unsigned int idkey); 

        //! insert Jet to the map. Once inserted, the JetMap own the Jet object
        Jet*   insert(Jet* jet);
        size_t erase(unsigned int idkey) {return _map.erase(idkey);}

  ConstIter begin()                  const {return _map.begin();}
  ConstIter find(unsigned int idkey) const {return _map.find(idkey);}
  ConstIter end()                    const {return _map.end();}

  Iter begin()                   {return _map.begin();}
  Iter  find(unsigned int idkey) {return _map.find(idkey);}
  Iter   end()                   {return _map.end();}
  
private:
  Jet::ALGO  _algo;                 //< algorithm used to reconstruct jets         
  float _par;                       //< algorithm parameter setting
  std::set<Jet::SRC> _src;      //< list of sources (clusters, towers, etc)
  typ_JetMap _map; //< jet algorithm output storage
    
  ClassDef(JetMap, 1);
};

#endif
