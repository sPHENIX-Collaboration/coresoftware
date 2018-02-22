#ifndef __TrkrHitSet_H__
#define __TrkrHitSet_H__

#include <TObject.h>
#include <g4main/PHG4Hit.h>
#include "TrkrDefUtil.h"

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif
#include <iostream>

/*
 * Virtual base class of hit container
 */
class TrkrHitSet : public TObject
{
 public:
  // ctor
  TrkrHitSet() {}
  // dtor
  virtual ~TrkrHitSet() {}
  // TObject functions
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Copy(TrkrHitSet& hit);
  virtual void Reset();
  virtual void print() const {};

  // Set/Get the key (ID) for this set of hits
  virtual void SetHitSetKey(const TrkrDefs::hitsetkey key) { return; }
  virtual TrkrDefs::hitsetkey GetHitSetKey() const { return 0; }
  // Set/Get the key (ID) for the truth mapping object assoc. to this hit set
  virtual void SetTruthMapKey(const uint64_t key) { return; }
  virtual uint64_t GetTruthMapKey() const { return ULLONG_MAX; }
  // Set/Get the key (ID) for the mapping between hits & clusters
  virtual void SetHitClusMapKey(const uint64_t key) { return; }
  virtual uint64_t GetHitClusMapKey() const { return TrkrDefs::CLUSKEYMAX; }
 private:
  ClassDef(TrkrHitSet, 1);
};

#endif