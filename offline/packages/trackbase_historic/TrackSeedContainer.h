#ifndef TRACKBASEHISTORIC_TRACKSEEDCONTAINER_H
#define TRACKBASEHISTORIC_TRACKSEEDCONTAINER_H

#include "TrackSeed.h"

#include <phool/PHObject.h>

#include <iostream>
#include <vector>


class TrackSeedContainer : public PHObject
{
 public:
  typedef std::vector<TrackSeed*> Container;
  typedef std::vector<TrackSeed*>::const_iterator ConstIter;
  typedef std::vector<TrackSeed*>::iterator Iter;
  
  ~TrackSeedContainer() override {}
  
  void identify(std::ostream& os = std::cout) const override
  {
    os << "TrackSeedContainer base class" << std::endl;
  }
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }
  
  virtual bool empty() const { return true; }
  virtual std::size_t size() const { return 0; }
  virtual void clear() {}

  virtual const TrackSeed* get(const std::size_t) const { return nullptr; }
  virtual TrackSeed* get(const std::size_t) { return nullptr; }
  virtual TrackSeed* insert(const TrackSeed*) { return nullptr; }
  virtual Iter erase(const std::size_t);

  virtual ConstIter begin() const;
  virtual ConstIter find(const std::size_t key) const;
  virtual ConstIter end() const;
  
  virtual Iter begin();
  virtual Iter find(const std::size_t key);
  virtual Iter end();

 protected:
  TrackSeedContainer() {}

 private:
  ClassDefOverride(TrackSeedContainer, 1);
};

#endif
