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
    int layer = TrkrDefs::getLayer(entry.first);
    os << "   hitset key: "  << entry.first << " layer " << layer
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

void 
TrkrHitTruthAssoc::findOrAddAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey, const PHG4HitDefs::keytype g4hitkey)
{
  // the association we want is between TrkrHit and PHG4Hit, but we need to know which TrkrHitSet the TrkrHit is in
  // check if this association already exists
  // We need all hitsets with this key

  std::pair<MMap::iterator, MMap::iterator> hitsetrange = m_map.equal_range(hitsetkey);
  MMap::iterator mapiter = hitsetrange.first;
  for (mapiter = hitsetrange.first; mapiter != hitsetrange.second; mapiter++)
    {
      if(mapiter->second.first == hitkey && mapiter->second.second == g4hitkey)
	{
	  // exists
	  //std::cout << "Association with hitsetkey " << hitsetkey << " hitkey " << hitkey << " g4hitkey " << g4hitkey << " already exists " << std::endl;
	  return;
	}
    }
  
  // Does not exist, create it
  std::pair<TrkrDefs::hitkey, PHG4HitDefs::keytype> assoc = std::make_pair(hitkey, g4hitkey);
  m_map.insert (std::make_pair(hitsetkey, assoc));
  //std::cout << "Added association with hitsetkey " << hitsetkey << " hitkey " << hitkey << " g4hitkey " << g4hitkey << std::endl;
}

TrkrHitTruthAssoc::ConstRange 
TrkrHitTruthAssoc::getCells(const TrkrDefs::hitsetkey hitsetkey, const unsigned int hidx)
{
return std::make_pair(m_map.end(), m_map.end());
}

