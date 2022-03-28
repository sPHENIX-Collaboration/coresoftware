#ifndef G4BBC_BBCVERTEXMAP_H
#define G4BBC_BBCVERTEXMAP_H

#include <phool/PHObject.h>

#include <cstddef>  // for size_t
#include <iostream>
#include <map>

class BbcVertex;

class BbcVertexMap : public PHObject
{
 public:
  typedef std::map<unsigned int, BbcVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, BbcVertex*>::iterator Iter;

  ~BbcVertexMap() override {}

  void identify(std::ostream& os = std::cout) const override { os << "BbcVertexMap base class" << std::endl; }
  int isValid() const override { return 0; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int /*idkey*/) const { return 0; }
  virtual void clear() {}

  virtual const BbcVertex* get(unsigned int /*idkey*/) const { return nullptr; }
  virtual BbcVertex* get(unsigned int /*idkey*/) { return nullptr; }
  virtual BbcVertex* insert(BbcVertex* /*vertex*/) { return nullptr; }
  virtual size_t erase(unsigned int /*idkey*/) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  BbcVertexMap() {}

 private:
  ClassDefOverride(BbcVertexMap, 1);
};

#endif  // G4BBC_BBCVERTEXMAP_H
