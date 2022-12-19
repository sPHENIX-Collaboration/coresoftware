/**
 * @file trackbase/RawHitSetContainer.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation for RawHitSetContainer
 */
#include "RawHitSetContainer.h"

#include <TSystem.h>

#include <cstdlib>
#include <iostream>


namespace
{
  RawHitSetContainer::Map dummy_map;
}

void RawHitSetContainer::Reset()
{
  std::cout << "RawHitSetContainer: Reset() not implemented by daughter class" << std::endl;
  gSystem->Exit(1);
}

RawHitSetContainer::ConstIterator
RawHitSetContainer::addHitSet(RawHitSet* /*newhit*/)
{ return dummy_map.cbegin(); }

RawHitSetContainer::ConstIterator
RawHitSetContainer::addHitSetSpecifyKey(const TrkrDefs::hitsetkey /*key*/, RawHitSet* /*newhit*/)
{ return dummy_map.cbegin(); }
 
RawHitSetContainer::Iterator
RawHitSetContainer::findOrAddHitSet(TrkrDefs::hitsetkey /*key*/)
{ return dummy_map.begin(); }

RawHitSetContainer::ConstRange
RawHitSetContainer::getHitSets(const TrkrDefs::TrkrId /*trackerid*/) const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }

RawHitSetContainer::ConstRange
RawHitSetContainer::getHitSets(const TrkrDefs::TrkrId, const uint8_t) const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }

RawHitSetContainer::ConstRange
RawHitSetContainer::getHitSets() const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }
