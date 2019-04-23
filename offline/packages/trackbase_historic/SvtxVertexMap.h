#ifndef __SVTXVERTEXMAP_H__
#define __SVTXVERTEXMAP_H__

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
  virtual SvtxVertexMap* Clone() const { return NULL; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int idkey) const { return 0; }
  virtual void clear() {}

  virtual const SvtxVertex* get(unsigned int idkey) const { return NULL; }
  virtual SvtxVertex* get(unsigned int idkey) { return NULL; }

  //! Add vertex to container. Note the container take to ownership
  virtual SvtxVertex* insert(const SvtxVertex* cluster) { return NULL; }
  //! legacy interface. Add vertex to container. Note the container does not take ownership
  virtual SvtxVertex* insert_clone(const SvtxVertex* vertex) { return NULL; }

  virtual size_t erase(unsigned int idkey) { return 0; }

  virtual ConstIter begin() const { return VertexMap().end(); }
  virtual ConstIter find(unsigned int idkey) const { return VertexMap().end(); }
  virtual ConstIter end() const { return VertexMap().end(); }

  virtual Iter begin() { return VertexMap().end(); }
  virtual Iter find(unsigned int idkey) { return VertexMap().end(); }
  virtual Iter end() { return VertexMap().end(); }

 protected:
  SvtxVertexMap() {}

 private:
  ClassDef(SvtxVertexMap, 1);
};

#endif
