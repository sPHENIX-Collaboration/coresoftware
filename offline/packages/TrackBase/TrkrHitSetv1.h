#ifndef TRACKBASE_TRKRHITSETV1_H
#define TRACKBASE_TRKRHITSETV1_H

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
  virtual void setHitSetKey(const TrkrDefs::hitsetkey key) { m_hitSetKey = key; }
  virtual TrkrDefs::hitsetkey getHitSetKey() const { return m_hitSetKey; }
  // Set/Get the key (ID) for the truth mapping object assoc. to this hit set
  virtual void setTruthMapKey(const uint64_t key) { m_truthMapKey = key; }
  virtual uint64_t getTruthMapKey() const { return m_truthMapKey; }
  // Set/Get the key (ID) for the mapping between hits & clusters
  virtual void setHitClusMapKey(const uint64_t key) { m_hitClusMapKey = key; }
  virtual uint64_t getHitClusMapKey() const { return m_hitClusMapKey; }
 protected:
 private:
  TrkrDefs::hitsetkey m_hitSetKey;
  uint64_t m_truthMapKey;
  uint64_t m_hitClusMapKey;
  ClassDef(TrkrHitSetv1,1);
};

#endif //TRACKBASE_TRKRHITSETV1_H
