#ifndef TRACKBASEHISTORIC_SVTXALIGNMENTSTATEMAP_H
#define TRACKBASEHISTORIC_SVTXALIGNMENTSTATEMAP_H

#include "SvtxAlignmentState.h"

#include <trackbase/TrkrDefs.h>

#include <iostream>
#include <map>

class SvtxAlignmentStateMap : public PHObject
{
 public:
  typedef std::vector<SvtxAlignmentState*> StateVec;
  typedef std::map<unsigned int, StateVec> StateMap;
  typedef StateMap::const_iterator ConstIter;
  typedef StateMap::iterator Iter;

  ~SvtxAlignmentStateMap() override {}

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxAlignmentStateMap base class " << std::endl;
  }
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  virtual bool empty() const { return true; }
  virtual std::size_t size() const { return 0; }
  virtual std::size_t count(unsigned int) const { return 0; }
  virtual void clear() {}

  virtual const StateVec get(unsigned int) const { return StateVec(); }
  virtual StateVec get(unsigned int) { return StateVec(); }
  virtual StateVec insertWithKey(unsigned int, StateVec) { return StateVec(); }
  virtual std::size_t erase(unsigned int) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(unsigned int) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(unsigned int);
  virtual Iter end();

 protected:
  SvtxAlignmentStateMap() {}

 private:
  ClassDefOverride(SvtxAlignmentStateMap, 1);
};

#endif
