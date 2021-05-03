#ifndef TRACKBASEHISTORIC_SVTXVERTEXMAP_H
#define TRACKBASEHISTORIC_SVTXVERTEXMAP_H

#include "SvtxVertex.h"

#include <phool/PHObject.h>
#include <iostream>
#include <map>

class SvtxVertexMap : public PHObject
{
 public:
  typedef std::map<unsigned int, SvtxVertex*> VertexMap;
  typedef std::map<unsigned int, SvtxVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, SvtxVertex*>::iterator Iter;

  virtual ~SvtxVertexMap() {}

  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "SvtxVertexMap base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual PHObject* CloneMe() const { return nullptr; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int idkey) const { return 0; }
  virtual void clear() {}

  virtual const SvtxVertex* get(unsigned int idkey) const { return nullptr; }
  virtual SvtxVertex* get(unsigned int idkey) { return nullptr; }

  //! Add vertex to container. Note the container take to ownership
  virtual SvtxVertex* insert(SvtxVertex* cluster) { return nullptr; }
  //! legacy interface. Add vertex to container. Note the container does not take ownership
  virtual SvtxVertex* insert_clone(const SvtxVertex* vertex) { return nullptr; }

  virtual size_t erase(unsigned int idkey) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  SvtxVertexMap() {}

 private:
  ClassDef(SvtxVertexMap, 1);
};

#endif
