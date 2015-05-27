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
   
  virtual ~JetMap() {}

  void identify(std::ostream &os = std::cout) const;
  void Reset() {}
  int  isValid() const {return 0;}

  // map content info ----------------------------------------------------------
  
  void      set_algo(Jet::ALGO algo) {return;}
  Jet::ALGO get_algo() const         {return Jet::NONE;}
  
  void  set_par(float par) {return;}
  float get_par() const    {return NAN;}

  // set access to list of source identifiers ----------------------------------
  
  bool empty_src() const        {return true;}
  void insert_src(Jet::SRC src) {return;}
 
  ConstSrcIter begin_src()             const {return std::set<Jet::SRC>().end();}
  ConstSrcIter  find_src(Jet::SRC src) const {return std::set<Jet::SRC>().end();}
  ConstSrcIter   end_src()             const {return std::set<Jet::SRC>().end();}

  SrcIter begin_src()             {return std::set<Jet::SRC>().end();}
  SrcIter  find_src(Jet::SRC src) {return std::set<Jet::SRC>().end();}
  SrcIter   end_src()             {return std::set<Jet::SRC>().end();}
  
  // map access to jets --------------------------------------------------------
 
  bool   empty()                   const {return true;}
  size_t  size()                   const {return 0;}
  size_t count(unsigned int idkey) const {return 0;}
  void   clear()                         {return;}
  
  const Jet* get(unsigned int idkey) const {return NULL;}
        Jet* get(unsigned int idkey)       {return NULL;} 

  Jet*   insert(Jet* jet)          {return NULL;}
  size_t erase(unsigned int idkey) {return 0;}

  ConstIter begin()                  const {return typ_JetMap().end();}
  ConstIter find(unsigned int idkey) const {return typ_JetMap().end();}
  ConstIter end()                    const {return typ_JetMap().end();}

  Iter begin()                   {return typ_JetMap().end();}
  Iter  find(unsigned int idkey) {return typ_JetMap().end();}
  Iter   end()                   {return typ_JetMap().end();}

protected:
  JetMap();
  
private:
    
  ClassDef(JetMap, 1);
};

#endif
