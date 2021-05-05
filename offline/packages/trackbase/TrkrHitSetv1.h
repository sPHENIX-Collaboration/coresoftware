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
#include <utility>           // for pair

// forward declaration
class TrkrHit;

class TrkrHitSetv1 : public TrkrHitSet
{
  
  public:
  
  TrkrHitSetv1() = default;
  
  virtual ~TrkrHitSetv1();

  virtual void identify(std::ostream& os = std::cout) const;
  
  virtual void Reset();

  virtual void setHitSetKey(const TrkrDefs::hitsetkey key) 
  { m_hitSetKey = key; }

  virtual TrkrDefs::hitsetkey getHitSetKey() const 
  { return m_hitSetKey; }

  virtual ConstIterator addHitSpecificKey(const TrkrDefs::hitkey, TrkrHit*);

  virtual void removeHit(TrkrDefs::hitkey);

  virtual TrkrHit* getHit(const TrkrDefs::hitkey) const;
  
  virtual ConstRange getHits() const;

  virtual unsigned int size() const
  { return m_hits.size(); }
   
  private:
  
  /// unique key for this object
  TrkrDefs::hitsetkey m_hitSetKey = TrkrDefs::HITSETKEYMAX;

  /// storage for TrkrHit objects
  Map m_hits;
  
  ClassDef(TrkrHitSetv1, 1);
};

#endif //TRACKBASE_TrkrHitSetv1_H
