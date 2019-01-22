/**
 * @file trackbase/TrkrHitTruthAssoc.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief Implementation of TrkrHitTruthAssoc
 */

#include "TrkrHitTruthAssoc.h"
#include "TrkrDefs.h"

TrkrHitTruthAssoc::TrkrHitTruthAssoc()
  : m_map()
{
}

TrkrHitTruthAssoc::~TrkrHitTruthAssoc()
{
}

void 
TrkrHitTruthAssoc::Reset()
{
m_map.clear();
}

void 
TrkrHitTruthAssoc::identify(std::ostream &os) const
{
  os << "-----TrkrHitTruthAssoc-----" << std::endl;
  os << "Number of associations: " << m_map.size() << std::endl;

  for ( auto& entry : m_map )
  {
    // os << "   cluster key: 0x" << std::hex << entry.first << std::dec
    os << "   hitset key: "  << entry.first 
       << " hit key: " << entry.second.first 
       << " g4hit key: " << entry.second.second 
       << std::endl;
  }

  os << "------------------------------" << std::endl;

  return;
}

void 
TrkrHitTruthAssoc::addAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey, const PHG4HitDefs::keytype g4hitkey)
{
// the association we want is between TrkrHit and PHG4Hit, but we need to know which TrkrHitSet the TrkrHit is in
std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> assoc = std::make_pair(hitkey, g4hitkey);
m_map.insert (std::make_pair(hitsetkey, assoc));
}

TrkrHitTruthAssoc::ConstRange 
TrkrHitTruthAssoc::getCells(const TrkrDefs::hitsetkey hitsetkey, const unsigned int hidx)
{
return std::make_pair(m_map.end(), m_map.end());
}

