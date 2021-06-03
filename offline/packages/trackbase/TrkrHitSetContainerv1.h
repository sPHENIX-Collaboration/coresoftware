#ifndef TRACKBASE_TrkrHitSetContainerv1_H
#define TRACKBASE_TrkrHitSetContainerv1_H
/**
 * @file trackbase/TrkrHitSetContainerv1.h
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 */

#include "TrkrDefs.h"       
#include "TrkrHitSetContainer.h"

#include <iostream>          // for cout, ostream
#include <map>
#include <utility>           // for pair

class TrkrHitSet;

/**
 * Container for TrkrHitSet objects
 */
class TrkrHitSetContainerv1 : public TrkrHitSetContainer
{
  
  public:

  TrkrHitSetContainerv1() = default;

  ~TrkrHitSetContainerv1() override
  { TrkrHitSetContainerv1::Reset(); }

  void Reset() override;

  void identify(std::ostream& = std::cout) const override;

  ConstIterator addHitSet(TrkrHitSet*) override;

  ConstIterator addHitSetSpecifyKey(const TrkrDefs::hitsetkey, TrkrHitSet*) override;

  void removeHitSet(TrkrDefs::hitsetkey ) override;

  void removeHitSet(TrkrHitSet* ) override;

  Iterator findOrAddHitSet(TrkrDefs::hitsetkey key) override;

  ConstRange getHitSets(const TrkrDefs::TrkrId trackerid) const override;

  ConstRange getHitSets(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const override;

  ConstRange getHitSets() const override;

  TrkrHitSet *findHitSet(TrkrDefs::hitsetkey key) override;

  unsigned int size() const override
  { return m_hitmap.size(); }

  private: 
  
  Map m_hitmap;
  
  ClassDefOverride(TrkrHitSetContainerv1, 1)
};

#endif //TRACKBASE_TrkrHitSetContainerv1_H
