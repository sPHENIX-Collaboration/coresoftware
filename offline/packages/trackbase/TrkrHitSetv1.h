#ifndef TRACKBASE_TRKRHITSETV1_H
#define TRACKBASE_TRKRHITSETV1_H

/**
 * @file trackbase/TrkrHitSetv1.h
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date 4 June 2018
 * @brief Container for storing TrkrHit's
 */
#include "TrkrDefs.h"
#include "TrkrHitSet.h"

#include <iostream>
#include <map>
#include <utility>  // for pair

// forward declaration
class TrkrHit;

class TrkrHitSetv1 : public TrkrHitSet
{
 public:
  TrkrHitSetv1() = default;

  ~TrkrHitSetv1() override
  {
    TrkrHitSetv1::Reset();
  }

  void identify(std::ostream& os = std::cout) const override;

  void Reset() override;

  void setHitSetKey(const TrkrDefs::hitsetkey key) override
  {
    m_hitSetKey = key;
  }

  TrkrDefs::hitsetkey getHitSetKey() const override
  {
    return m_hitSetKey;
  }

  ConstIterator addHitSpecificKey(const TrkrDefs::hitkey, TrkrHit*) override;

  void removeHit(TrkrDefs::hitkey) override;

  TrkrHit* getHit(const TrkrDefs::hitkey) const override;

  ConstRange getHits() const override;

  unsigned int size() const override
  {
    return m_hits.size();
  }

 private:
  /// unique key for this object
  TrkrDefs::hitsetkey m_hitSetKey = TrkrDefs::HITSETKEYMAX;

  /// storage for TrkrHit objects
  Map m_hits;

  ClassDefOverride(TrkrHitSetv1, 1);
};

#endif  //TRACKBASE_TrkrHitSetv1_H
