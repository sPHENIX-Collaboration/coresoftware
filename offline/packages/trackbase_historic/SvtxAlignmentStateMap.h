#ifndef TRACKBASEHISTORIC_SVTXALIGNMENTSTATEMAP_H
#define TRACKBASEHISTORIC_SVTXALIGNMENTSTATEMAP_H

#include "SvtxAlignmentState.h"

#include <trackbase/TrkrDefs.h>

#include <map>
#include <iostream>

class SvtxAlignmentStateMap : public PHObject
{
 public:
  typedef std::map<TrkrDefs::cluskey, SvtxAlignmentState*> StateMap;
  typedef StateMap::const_iterator ConstIter;
  typedef StateMap::iterator Iter;

  ~SvtxAlignmentStateMap() override {}

  void identify(std::ostream& os = std::cout ) const override
  {
    os << "SvtxAlignmentStateMap base class " << std::endl;
  }
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  virtual bool empty() const { return true; }
  virtual std::size_t size() const { return 0; }
  virtual std::size_t count(TrkrDefs::cluskey) const { return 0; }
  virtual void clear() {}

  virtual const SvtxAlignmentState* get(TrkrDefs::cluskey) const { return nullptr; }
  virtual SvtxAlignmentState* get(TrkrDefs::cluskey) { return nullptr; }
  virtual SvtxAlignmentState* insertWithKey(TrkrDefs::cluskey, SvtxAlignmentState*) { return nullptr; }
  virtual std::size_t erase(TrkrDefs::cluskey) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(TrkrDefs::cluskey cluskey) const;
  virtual ConstIter end() const;
  
  virtual Iter begin();
  virtual Iter find(TrkrDefs::cluskey cluskey);
  virtual Iter end();
  
 protected:
  SvtxAlignmentStateMap() {}

 private:
  ClassDefOverride(SvtxAlignmentStateMap, 1);

};

#endif
