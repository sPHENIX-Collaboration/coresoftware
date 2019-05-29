/**
 * @file trackbase/TrkrHitCellAssoc.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrHitCellAssoc
 */

#include "TrkrHitCellAssoc.h"

#include <g4detectors/PHG4CellDefs.h>  // for keytype

#include <iosfwd>                      // for ostream
#include <type_traits>                 // for __decay_and_strip<>::__type

TrkrHitCellAssoc::TrkrHitCellAssoc()
  : m_map()
{
}

TrkrHitCellAssoc::~TrkrHitCellAssoc()
{
}

void 
TrkrHitCellAssoc::Reset()
{
  m_map.clear();
}

void 
TrkrHitCellAssoc::identify(std::ostream &os) const
{
}

void 
TrkrHitCellAssoc::addAssoc(const TrkrDefs::hitsetkey hset, const unsigned int hidx, const PHG4CellDefs::keytype ckey)
{
}

TrkrHitCellAssoc::ConstRange 
TrkrHitCellAssoc::getCells(const TrkrDefs::hitsetkey hset, const unsigned int hidx)
{
  return std::make_pair(m_map.end(), m_map.end());
}

