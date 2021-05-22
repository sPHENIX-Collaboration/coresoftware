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

  virtual ~TrkrHitSetv1()
  {
    TrkrHitSetv1::Reset();
  }

  virtual void identify(std::ostream& os = std::cout) const override;

  virtual void Reset() override;

  virtual void setHitSetKey(const TrkrDefs::hitsetkey key) override
  {
    m_hitSetKey = key;
  }

  virtual TrkrDefs::hitsetkey getHitSetKey() const override
  {
    return m_hitSetKey;
  }

  virtual ConstIterator addHitSpecificKey(const TrkrDefs::hitkey, TrkrHit*) override;

  virtual void removeHit(TrkrDefs::hitkey) override;

  virtual TrkrHit* getHit(const TrkrDefs::hitkey) const override;

  virtual ConstRange getHits() const override;

  virtual unsigned int size() const override
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
