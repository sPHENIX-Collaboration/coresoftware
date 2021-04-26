// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4VERTEX_GLOBALVERTEXMAP_H
#define G4VERTEX_GLOBALVERTEXMAP_H

#include <phool/PHObject.h>

#include <iostream>
#include <map>

class GlobalVertex;

class GlobalVertexMap : public PHObject
{
 public:
  typedef std::map<unsigned int, GlobalVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, GlobalVertex*>::iterator Iter;

  virtual ~GlobalVertexMap() {}

  virtual void identify(std::ostream& os = std::cout) const { os << "GlobalVertexMap base class" << std::endl; }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int idkey) const { return 0; }
  virtual void clear() {}

  virtual const GlobalVertex* get(unsigned int idkey) const { return nullptr; }
  virtual GlobalVertex* get(unsigned int idkey) { return nullptr; }
  virtual GlobalVertex* insert(GlobalVertex* vertex) { return nullptr; }
  virtual size_t erase(unsigned int idkey) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  GlobalVertexMap() {}

 private:
  ClassDef(GlobalVertexMap, 1);
};

#endif  // G4VERTEX_GLOBALVERTEXMAP_H
