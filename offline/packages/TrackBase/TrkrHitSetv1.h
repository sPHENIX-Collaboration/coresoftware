#ifndef __TrkrHitSetV1_H__
#define __TrkrHitSetV1_H__

#include "TrkrDefUtil.h"
#include "TrkrHitSet.h"

#ifdef __CINT__
#include <stdint.h>
#else
#include <cstdint>
#endif
#include <iostream>
#include <map>

/**
 * Version 1 of hit container
 */
class TrkrHitSetv1 : public TrkrHitSet
{
 public:
  // ctor
  TrkrHitSetv1();

  // dtor
  virtual ~TrkrHitSetv1() {}
  // TObject classes
  virtual void identify(std::ostream& os = std::cout) const;
  virtual void Reset();
  virtual void print() const;

  // Set/Get the key (ID) for this set of hits
  virtual void SetHitSetKey(const TrkrDefs::hitsetkey key) { hitset_key_ = key; }
  virtual TrkrDefs::hitsetkey GetHitSetKey() const { return hitset_key_; }
  // Set/Get the key (ID) for the truth mapping object assoc. to this hit set
  virtual void SetTruthMapKey(const uint64_t key) { truth_map_key_ = key; }
  virtual uint64_t GetTruthMapKey() const { return truth_map_key_; }
  // Set/Get the key (ID) for the mapping between hits & clusters
  virtual void SetHitClusMapKey(const uint64_t key) { hit_clus_map_key_ = key; }
  virtual uint64_t GetHitClusMapKey() const { return hit_clus_map_key_; }
 protected:
 private:
  TrkrDefs::hitsetkey hitset_key_;
  uint64_t truth_map_key_;
  uint64_t hit_clus_map_key_;
  ClassDef(TrkrHitSetv1,1);
};

#endif