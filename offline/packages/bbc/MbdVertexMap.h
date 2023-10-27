#ifndef G4MBD_MBDVERTEXMAP_H
#define G4MBD_MBDVERTEXMAP_H

#include <phool/PHObject.h>

#include <cstddef>  // for size_t
#include <iostream>
#include <map>

class MbdVertex;

class MbdVertexMap : public PHObject
{
 public:
  typedef std::map<unsigned int, MbdVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, MbdVertex*>::iterator Iter;

  ~MbdVertexMap() override {}

  void identify(std::ostream& os = std::cout) const override { os << "MbdVertexMap base class" << std::endl; }
  int isValid() const override { return 0; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int /*idkey*/) const { return 0; }
  virtual void clear() {}

  virtual const MbdVertex* get(unsigned int /*idkey*/) const { return nullptr; }
  virtual MbdVertex* get(unsigned int /*idkey*/) { return nullptr; }
  virtual MbdVertex* insert(MbdVertex* /*vertex*/) { return nullptr; }
  virtual size_t erase(unsigned int /*idkey*/) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  MbdVertexMap() {}

 private:
  ClassDefOverride(MbdVertexMap, 1);
};

#endif  // G4MBD_MBDVERTEXMAP_H
