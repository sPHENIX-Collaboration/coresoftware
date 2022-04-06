#ifndef TRACKBASEHISTORIC_SVTXTRACKSEEDCONTAINER_H
#define TRACKBASEHISTORIC_SVTXTRACKSEEDCONTAINER_H

#include "SvtxTrackSeed.h"

#include <phool/PHObject.h>

#include <iostream>
#include <vector>


class SvtxTrackSeedContainer : public PHObject
{
 public:
  typedef std::vector<SvtxTrackSeed*> TrackSeedContainer;
  typedef std::vector<SvtxTrackSeed*>::const_iterator ConstIter;
  typedef std::vector<SvtxTrackSeed*>::iterator Iter;
  
  ~SvtxTrackSeedContainer() override {}
  
  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxTrackSeedContainer base class" << std::endl;
  }
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }
  
  virtual bool empty() const { return true; }
  virtual std::size_t size() const { return 0; }
  virtual void clear() {}

  virtual const SvtxTrackSeed* get(const std::size_t) const { return nullptr; }
  virtual SvtxTrackSeed* get(const std::size_t) { return nullptr; }
  virtual SvtxTrackSeed* insert(const SvtxTrackSeed*) { return nullptr; }
  virtual Iter erase(const std::size_t);

  virtual ConstIter begin() const;
  virtual ConstIter find(const std::size_t key) const;
  virtual ConstIter end() const;
  
  virtual Iter begin();
  virtual Iter find(const std::size_t key);
  virtual Iter end();

 protected:
  SvtxTrackSeedContainer() {}

 private:
  ClassDefOverride(SvtxTrackSeedContainer, 1);
};

#endif
