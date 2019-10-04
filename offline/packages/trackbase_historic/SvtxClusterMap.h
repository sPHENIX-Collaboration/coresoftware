#ifndef __SVTXCLUSTERMAP_H__
#define __SVTXCLUSTERMAP_H__

#include "SvtxCluster.h"

#include <phool/PHObject.h>
#include <iostream>
#include <map>

class SvtxClusterMap : public PHObject
{
 public:
  typedef std::map<unsigned int, SvtxCluster*> ClusterMap;
  typedef std::map<unsigned int, SvtxCluster*>::const_iterator ConstIter;
  typedef std::map<unsigned int, SvtxCluster*>::iterator Iter;

  virtual ~SvtxClusterMap() {}

  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "SvtxClusterMap base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual PHObject* CloneMe() const { return nullptr; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int idkey) const { return 0; }
  virtual void clear() {}

  virtual const SvtxCluster* get(unsigned int idkey) const { return nullptr; }
  virtual SvtxCluster* get(unsigned int idkey) { return nullptr; }
  virtual SvtxCluster* insert(const SvtxCluster* cluster) { return nullptr; }
  virtual size_t erase(unsigned int idkey) { return 0; }

  virtual ConstIter begin() const { return ClusterMap().end(); }
  virtual ConstIter find(unsigned int idkey) const { return ClusterMap().end(); }
  virtual ConstIter end() const { return ClusterMap().end(); }

  virtual Iter begin() { return ClusterMap().end(); }
  virtual Iter find(unsigned int idkey) { return ClusterMap().end(); }
  virtual Iter end() { return ClusterMap().end(); }

 protected:
  SvtxClusterMap() {}

 private:
  ClassDef(SvtxClusterMap, 1);
};

#endif
