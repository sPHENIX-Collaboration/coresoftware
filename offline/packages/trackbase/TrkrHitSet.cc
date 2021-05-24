/**
 * @file trackbase/TrkrHitSet.cc
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * @brief Implementation of TrkrHitSet
 */
#include "TrkrHitSet.h"

namespace
{
  TrkrHitSet::Map dummy_map;
}

TrkrHitSet::ConstIterator
TrkrHitSet::addHitSpecificKey(const TrkrDefs::hitkey, TrkrHit*)
{ return dummy_map.cbegin(); }

TrkrHitSet::ConstRange
TrkrHitSet::getHits() const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }
