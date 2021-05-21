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

  virtual ~TrkrHitSetContainerv1()
  { TrkrHitSetContainerv1::Reset(); }

  virtual void Reset() override;

  virtual void identify(std::ostream& = std::cout) const override;

  virtual ConstIterator addHitSet(TrkrHitSet*) override;

  virtual ConstIterator addHitSetSpecifyKey(const TrkrDefs::hitsetkey, TrkrHitSet*) override;

  virtual void removeHitSet(TrkrDefs::hitsetkey ) override;

  virtual void removeHitSet(TrkrHitSet* ) override;

  virtual Iterator findOrAddHitSet(TrkrDefs::hitsetkey key) override;

  virtual ConstRange getHitSets(const TrkrDefs::TrkrId trackerid) const override;

  virtual ConstRange getHitSets(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const override;

  virtual ConstRange getHitSets() const override;

  virtual TrkrHitSet *findHitSet(TrkrDefs::hitsetkey key) override;

  virtual unsigned int size() const override
  { return m_hitmap.size(); }

  private: 
  
  Map m_hitmap;
  
  ClassDefOverride(TrkrHitSetContainerv1, 1)
};

#endif //TRACKBASE_TrkrHitSetContainerv1_H
