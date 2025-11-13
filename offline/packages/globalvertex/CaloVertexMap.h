// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GLOBALVERTEX_CALOVERTEXMAP_H
#define GLOBALVERTEX_CALOVERTEXMAP_H

#include <phool/PHObject.h>

#include <cstddef>  // for size_t
#include <iostream>
#include <map>

class CaloVertex;

class CaloVertexMap : public PHObject
{
 public:
  typedef std::map<unsigned int, CaloVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, CaloVertex*>::iterator Iter;

  ~CaloVertexMap() override = default;

  void identify(std::ostream& os = std::cout) const override { os << "CaloVertexMap base class" << std::endl; }
  int isValid() const override { return 0; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int /*idkey*/) const { return 0; }
  virtual void clear() {}

  virtual const CaloVertex* get(unsigned int /*idkey*/) const { return nullptr; }
  virtual CaloVertex* get(unsigned int /*idkey*/) { return nullptr; }
  virtual CaloVertex* insert(CaloVertex* /*vertex*/) { return nullptr; }
  virtual size_t erase(unsigned int /*idkey*/) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  CaloVertexMap() = default;

 private:
  ClassDefOverride(CaloVertexMap, 1);
};

#endif  // G4CALO_CALOVERTEXMAP_H
