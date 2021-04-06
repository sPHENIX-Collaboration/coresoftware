#ifndef TRACKBASEHISTORIC_SVTXTRACKMAP_H
#define TRACKBASEHISTORIC_SVTXTRACKMAP_H

#include "SvtxTrack.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>

class SvtxTrackMap : public PHObject
{
 public:
  typedef std::map<unsigned int, SvtxTrack*> TrackMap;
  typedef std::map<unsigned int, SvtxTrack*>::const_iterator ConstIter;
  typedef std::map<unsigned int, SvtxTrack*>::iterator Iter;

  virtual ~SvtxTrackMap() {}

  virtual void identify(std::ostream& os = std::cout) const
  {
    os << "SvtxTrackMap base class" << std::endl;
  }
  virtual void Reset() {}
  virtual int isValid() const { return 0; }
  virtual PHObject* CloneMe() const { return nullptr; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual size_t count(unsigned int idkey) const { return 0; }
  virtual void clear() {}

  virtual const SvtxTrack* get(unsigned int idkey) const { return nullptr; }
  virtual SvtxTrack* get(unsigned int idkey) { return nullptr; }
  virtual SvtxTrack* insert(const SvtxTrack* cluster) { return nullptr; }
  virtual size_t erase(unsigned int idkey) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int idkey) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int idkey);
  virtual Iter end();

 protected:
  SvtxTrackMap() {}

 private:
  ClassDef(SvtxTrackMap, 1);
};

#endif
