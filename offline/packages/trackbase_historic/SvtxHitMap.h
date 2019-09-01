#ifndef __SVTXHITMAP_H__
#define __SVTXHITMAP_H__

#include "SvtxHit.h"

#include <phool/PHObject.h>
#include <iostream>
#include <map>

class SvtxHitMap : public PHObject
{
 public:
  typedef std::map<unsigned int, SvtxHit*> HitMap;
  typedef std::map<unsigned int, SvtxHit*>::const_iterator ConstIter;
  typedef std::map<unsigned int, SvtxHit*>::iterator Iter;

  virtual ~SvtxHitMap() {}

  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "SvtxHitMap base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual PHObject* CloneMe() const { return nullptr; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int idkey) const { return 0; }
  virtual void clear() {}

  virtual const SvtxHit* get(unsigned int idkey) const { return nullptr; }
  virtual SvtxHit* get(unsigned int idkey) { return nullptr; }
  virtual SvtxHit* insert(const SvtxHit* hit) { return nullptr; }
  virtual size_t erase(unsigned int idkey) { return 0; }

  virtual ConstIter begin() const { return HitMap().end(); }
  virtual ConstIter find(unsigned int idkey) const { return HitMap().end(); }
  virtual ConstIter end() const { return HitMap().end(); }

  virtual Iter begin() { return HitMap().end(); }
  virtual Iter find(unsigned int idkey) { return HitMap().end(); }
  virtual Iter end() { return HitMap().end(); }

 protected:
  SvtxHitMap() {}

 private:
  ClassDef(SvtxHitMap, 1);
};

#endif
