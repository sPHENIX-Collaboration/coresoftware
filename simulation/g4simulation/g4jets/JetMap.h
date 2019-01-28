#ifndef G4JET_JETMAP_H
#define G4JET_JETMAP_H

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

  virtual void    identify(std::ostream &os = std::cout) const;
  virtual void    Reset() {}
  virtual int     isValid() const {return 0;}
  virtual JetMap* Clone() const {return nullptr;}

  // map content info ----------------------------------------------------------
  
  virtual void      set_algo(Jet::ALGO algo) {return;}
  virtual Jet::ALGO get_algo() const         {return Jet::NONE;}
  
  virtual void  set_par(float par) {return;}
  virtual float get_par() const    {return NAN;}

  // set access to list of source identifiers ----------------------------------
  
  virtual bool empty_src() const        {return true;}
  virtual void insert_src(Jet::SRC src) {return;}
 
  virtual ConstSrcIter begin_src()             const {return std::set<Jet::SRC>().end();}
  virtual ConstSrcIter  find_src(Jet::SRC src) const {return std::set<Jet::SRC>().end();}
  virtual ConstSrcIter   end_src()             const {return std::set<Jet::SRC>().end();}

  virtual SrcIter begin_src()             {return std::set<Jet::SRC>().end();}
  virtual SrcIter  find_src(Jet::SRC src) {return std::set<Jet::SRC>().end();}
  virtual SrcIter   end_src()             {return std::set<Jet::SRC>().end();}
  
  // map access to jets --------------------------------------------------------
 
  virtual bool   empty()                   const {return true;}
  virtual size_t  size()                   const {return 0;}
  virtual size_t count(unsigned int idkey) const {return 0;}
  virtual void   clear()                         {return;}
  
  virtual const Jet* get(unsigned int idkey) const {return nullptr;}
  virtual       Jet* get(unsigned int idkey)       {return nullptr;} 

  virtual Jet*   insert(Jet* jet)          {return nullptr;}
  virtual size_t erase(unsigned int idkey) {return 0;}

  virtual ConstIter begin()                  const {return typ_JetMap().end();}
  virtual ConstIter find(unsigned int idkey) const {return typ_JetMap().end();}
  virtual ConstIter end()                    const {return typ_JetMap().end();}

  virtual Iter begin()                   {return typ_JetMap().end();}
  virtual Iter  find(unsigned int idkey) {return typ_JetMap().end();}
  virtual Iter   end()                   {return typ_JetMap().end();}

protected:
  JetMap();
  
private:
    
  ClassDef(JetMap, 1);
};

#endif
