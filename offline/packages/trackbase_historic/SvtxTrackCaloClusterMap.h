#ifndef TRACKBASEHISTORIC_SVTXTRACKCALOCLUSTERMAP_H
#define TRACKBASEHISTORIC_SVTXTRACKCALOCLUSTERMAP_H

#include <calobase/RawCluster.h>
#include "SvtxTrack.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>

class SvtxTrackCaloClusterMap : public PHObject
{
 public:
  typedef std::map<SvtxTrack*, std::vector<RawCluster*>> Map;
  typedef std::map<SvtxTrack*, std::vector<RawCluster*>>::const_iterator ConstIter;
  typedef std::map<SvtxTrack*, std::vector<RawCluster*>>::iterator Iter;

  ~SvtxTrackCaloClusterMap() override {}

  void identify(std::ostream& os = std::cout) const override
  {
    os << "SvtxTrackCaloClusterMap base class" << std::endl;
  }
  int isValid() const override { return 0; }
  PHObject* CloneMe() const override { return nullptr; }

  virtual bool empty() const { return true; }
  virtual size_t size() const { return 0; }
  virtual void clear() {}

  virtual const std::vector<RawCluster*> get(SvtxTrack*) const
  {
    std::vector<RawCluster*> dum;
    return dum;
  }
  virtual std::vector<RawCluster*> get(SvtxTrack*)
  {
    std::vector<RawCluster*> dum;
    return dum;
  }
  virtual std::vector<RawCluster*> insert(SvtxTrack*, std::vector<RawCluster*>)
  {
    std::vector<RawCluster*> dummy;
    return dummy;
  }
  virtual std::vector<RawCluster*> insert(SvtxTrack*, RawCluster*)
  {
    std::vector<RawCluster*> dummy;
    return dummy;
  }
  virtual size_t erase(SvtxTrack*) { return 0; }

  virtual ConstIter begin() const;
  virtual ConstIter find(SvtxTrack*) const;
  virtual ConstIter end() const;

  virtual Iter begin();
  virtual Iter find(SvtxTrack*);
  virtual Iter end();

 protected:
  SvtxTrackCaloClusterMap() {}

 private:
  ClassDefOverride(SvtxTrackCaloClusterMap, 1);
};

#endif
