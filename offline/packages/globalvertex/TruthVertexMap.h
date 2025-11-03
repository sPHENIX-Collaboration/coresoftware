// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_TRUTHVERTEXMAP_H
#define GLOBALVERTEX_TRUTHVERTEXMAP_H

#include <phool/PHObject.h>

#include <cstddef>  // for size_t
#include <iostream>
#include <map>

class TruthVertex;

class TruthVertexMap : public PHObject
{
 public:
  typedef std::map<unsigned int, TruthVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, TruthVertex*>::iterator Iter;

  ~TruthVertexMap() override = default;

  void identify(std::ostream& os = std::cout) const override { os << "TruthVertexMap base class" << std::endl; }
  int isValid() const override { return 0; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int /*idkey*/) const { return 0; }
  virtual void clear() {}

  virtual const TruthVertex* get(unsigned int /*idkey*/) const { return nullptr; }
  virtual TruthVertex* get(unsigned int /*idkey*/) { return nullptr; }
  virtual TruthVertex* insert(TruthVertex* /*vertex*/) { return nullptr; }
  virtual size_t erase(unsigned int /*idkey*/) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  TruthVertexMap() = default;

 private:
  ClassDefOverride(TruthVertexMap, 1);
};

#endif  // G4TRUTH_TRUTHVERTEXMAP_H