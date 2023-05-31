#ifndef TRACKBASE_TrkrHitSetContainerv2_H
#define TRACKBASE_TrkrHitSetContainerv2_H

#include "TrkrDefs.h"
#include "TrkrHitSetContainer.h"

#include <TClonesArray.h>
#include <iostream>  // for cout, ostream
#include <map>
#include <string>   // for pair
#include <utility>  // for pair

class TrkrHitSet;

/**
 * IO and memory efficient container using TClonesArray
 * TrkrHitSet used in this container must impliment TrkrHitSet::Clear() function for fast reset without calling ~TrkrHitSet()
 */
class TrkrHitSetContainerv2 final : public TrkrHitSetContainer
{
public:
  //! only used in ROOT IO. Do NOT use this constructor in user code
  TrkrHitSetContainerv2() = default;

  TrkrHitSetContainerv2(const std::string& hitsetclass, const size_t estimated_size);

  ~TrkrHitSetContainerv2() override
  {
  }

  void Reset() override;

  void identify(std::ostream& = std::cout) const override;

  ConstIterator addHitSet(TrkrHitSet*) override;

  ConstIterator addHitSetSpecifyKey(const TrkrDefs::hitsetkey, TrkrHitSet*) override;

  //! Add a TrkrHitSet to TrkrHitSetContainerv2, return the iterater to the newly constructed TrkrHitSet
  ConstIterator addHitSetSpecifyKey(const TrkrDefs::hitsetkey key);

  void removeHitSet(TrkrDefs::hitsetkey) override;

  void removeHitSet(TrkrHitSet*) override;

  Iterator findOrAddHitSet(TrkrDefs::hitsetkey key) override;

  ConstRange getHitSets(const TrkrDefs::TrkrId trackerid) const override;

  ConstRange getHitSets(const TrkrDefs::TrkrId trackerid, const uint8_t layer) const override;

  ConstRange getHitSets() const override;

  TrkrHitSet* findHitSet(TrkrDefs::hitsetkey key) override;

  unsigned int size() const override
  {
    return m_hitArray.GetEntriesFast();
  }

 private:
  //! endure m_hitmap and m_hitArray are in sync, in particular m_hitmap need rebuild after reading back from the DST
  void syncMapArray(void) const;

  //! used for indexing only, not used in storage. syncMapArray(void) rebuilds the map after DST readback
  mutable Map m_hitmap;  //!

  TClonesArray m_hitArray;

  ClassDefOverride(TrkrHitSetContainerv2, 1)
};

#endif  // TRACKBASE_TrkrHitSetContainerv2_H
