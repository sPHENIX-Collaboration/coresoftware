#ifndef TRACKBASE_TRKRHITSET_H
#define TRACKBASE_TRKRHITSET_H

#include "TrkrDefUtil.h"

#include <g4main/PHG4Hit.h>

#include <TObject.h>

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
  virtual void setHitSetKey(const TrkrDefs::hitsetkey key) { return; }
  virtual TrkrDefs::hitsetkey getHitSetKey() const { return 0; }
  // Set/Get the key (ID) for the truth mapping object assoc. to this hit set
  virtual void setTruthMapKey(const uint64_t key) { return; }
  virtual uint64_t getTruthMapKey() const { return ULLONG_MAX; }
  // Set/Get the key (ID) for the mapping between hits & clusters
  virtual void setHitClusMapKey(const uint64_t key) { return; }
  virtual uint64_t getHitClusMapKey() const { return TrkrDefs::CLUSKEYMAX; }
 private:
  ClassDef(TrkrHitSet, 1);
};

#endif //TRACKBASE_TRKRHITSET_H
