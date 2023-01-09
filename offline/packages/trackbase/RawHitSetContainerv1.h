#ifndef TRACKBASE_RawHitSetContainerv1_H
#define TRACKBASE_RawHitSetContainerv1_H
/**
 * @file trackbase/RawHitSetContainerv1.h
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 */

#include "TrkrDefs.h"       
#include "RawHitSetContainer.h"

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair

class RawHitSet;

/**
 * Container for RawHitSet objects
 */
class RawHitSetContainerv1 : public RawHitSetContainer
{
  
  public:

  RawHitSetContainerv1() = default;

  ~RawHitSetContainerv1() override
  { RawHitSetContainerv1::Reset(); }

  void Reset() override;

  void identify(std::ostream& = std::cout) const override;

  ConstIterator addHitSet(RawHitSet*) override;

  ConstIterator addHitSetSpecifyKey(const TrkrDefs::hitsetkey, RawHitSet*) override;

  void removeHitSet(TrkrDefs::hitsetkey ) override;

  void removeHitSet(RawHitSet* ) override;

  Iterator findOrAddHitSet(TrkrDefs::hitsetkey key) override;

  ConstRange getHitSets(const TrkrDefs::TrkrId trackerid) const override;

  ConstRange getHitSets(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const override;

  ConstRange getHitSets() const override;

  RawHitSet *findHitSet(TrkrDefs::hitsetkey key) override;

  unsigned int size() const override
  { return m_hitmap.size(); }

  private: 
  
  Map m_hitmap;
  
  ClassDefOverride(RawHitSetContainerv1, 1)
};

#endif //TRACKBASE_RawHitSetContainerv1_H
