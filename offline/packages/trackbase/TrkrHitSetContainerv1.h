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
  { Reset(); }

  virtual void Reset();

  virtual void identify(std::ostream& = std::cout) const;

  virtual ConstIterator addHitSet(TrkrHitSet*);

  virtual ConstIterator addHitSetSpecifyKey(const TrkrDefs::hitsetkey, TrkrHitSet*);

  virtual void removeHitSet(TrkrDefs::hitsetkey );

  virtual void removeHitSet(TrkrHitSet* );

  virtual Iterator findOrAddHitSet(TrkrDefs::hitsetkey key);

  virtual ConstRange getHitSets(const TrkrDefs::TrkrId trackerid) const;

  virtual ConstRange getHitSets(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const;

  virtual ConstRange getHitSets() const;

  virtual TrkrHitSet *findHitSet(TrkrDefs::hitsetkey key);

  virtual unsigned int size() const
  { return m_hitmap.size(); }

  private: 
  
  Map m_hitmap;
  
  ClassDef(TrkrHitSetContainerv1, 1)
};

#endif //TRACKBASE_TrkrHitSetContainerv1_H
