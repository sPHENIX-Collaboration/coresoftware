/**
 * @file trackbase/TrkrHitSetContainer.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation for TrkrHitSetContainer
 */
#include "TrkrHitSetContainer.h"

#include <TSystem.h>

#include <cstdlib>
#include <iostream>


namespace
{
  TrkrHitSetContainer::Map dummy_map;
}

void TrkrHitSetContainer::Reset()
{
  std::cout << "TrkrHitSetContainer: Reset() not implemented by daughter class" << std::endl;
  gSystem->Exit(1);
}

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::addHitSet(TrkrHitSet* /*newhit*/)
{ return dummy_map.cbegin(); }

TrkrHitSetContainer::ConstIterator
TrkrHitSetContainer::addHitSetSpecifyKey(const TrkrDefs::hitsetkey /*key*/, TrkrHitSet* /*newhit*/)
{ return dummy_map.cbegin(); }
 
TrkrHitSetContainer::Iterator
TrkrHitSetContainer::findOrAddHitSet(TrkrDefs::hitsetkey /*key*/)
{ return dummy_map.begin(); }

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::getHitSets(const TrkrDefs::TrkrId /*trackerid*/) const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::getHitSets(const TrkrDefs::TrkrId, const uint8_t) const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }

TrkrHitSetContainer::ConstRange
TrkrHitSetContainer::getHitSets() const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }
