#ifndef G4JET_JETMAP_H
#define G4JET_JETMAP_H

#include "Jet.h"

#include <phool/PHObject.h>

#include <cmath>
#include <cstddef>          // for size_t
#include <iostream>
#include <map>
#include <set>

class JetMap : public PHObject
{
 public:
  // jet object iterators
  typedef std::map<unsigned int, Jet*> typ_JetMap;
  typedef typ_JetMap::const_iterator ConstIter;
  typedef typ_JetMap::iterator Iter;

  // source identifier iterators
  typedef std::set<Jet::SRC>::const_iterator ConstSrcIter;
  typedef std::set<Jet::SRC>::iterator SrcIter;

  JetMap(){}
  virtual ~JetMap() {}

  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual PHObject* CloneMe() const { return nullptr; }

  // map content info ----------------------------------------------------------

  virtual void set_algo(Jet::ALGO algo) { return; }
  virtual Jet::ALGO get_algo() const { return Jet::NONE; }

  virtual void set_par(float par) { return; }
  virtual float get_par() const { return NAN; }

  // set access to list of source identifiers ----------------------------------

  virtual bool empty_src() const { return true; }
  virtual void insert_src(Jet::SRC src) { return; }

  virtual ConstSrcIter begin_src() const;
  virtual ConstSrcIter find_src(Jet::SRC src) const;
  virtual ConstSrcIter end_src() const;

  virtual SrcIter begin_src();
  virtual SrcIter find_src(Jet::SRC src);
  virtual SrcIter end_src();

  // map access to jets --------------------------------------------------------

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int idkey) const { return 0; }
  virtual void clear() { return; }

  virtual const Jet* get(unsigned int idkey) const { return nullptr; }
  virtual Jet* get(unsigned int idkey) { return nullptr; }

  virtual Jet* insert(Jet* jet) { return nullptr; }
  virtual size_t erase(unsigned int idkey) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 private:
  ClassDef(JetMap, 1);
};

#endif
