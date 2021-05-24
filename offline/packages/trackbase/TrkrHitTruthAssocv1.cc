/**
 * @file trackbase/TrkrHitTruthAssocv1.cc
 * @author D. McGlinchey, H. PEREIRA DA COSTA
 * @date June 2018
 * @brief Implementation of TrkrHitTruthAssocv1
 */

#include "TrkrHitTruthAssocv1.h"
#include "TrkrDefs.h"

#include <g4main/PHG4HitDefs.h>

#include <algorithm>
#include <ostream>               // for operator<<, endl, basic_ostream, bas...

void
TrkrHitTruthAssocv1::Reset()
{ m_map.clear(); }

void
TrkrHitTruthAssocv1::identify(std::ostream &os) const
{
  os << "-----TrkrHitTruthAssocv1-----" << std::endl;
  os << "Number of associations: " << m_map.size() << std::endl;

  for( const auto& entry : m_map )
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
TrkrHitTruthAssocv1::addAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey, const PHG4HitDefs::keytype g4hitkey)
{
  // the association we want is between TrkrHit and PHG4Hit, but we need to know which TrkrHitSet the TrkrHit is in
  m_map.insert (std::make_pair(hitsetkey, std::make_pair(hitkey, g4hitkey)));
}

void
TrkrHitTruthAssocv1::findOrAddAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey, const PHG4HitDefs::keytype g4hitkey)
{
  // the association we want is between TrkrHit and PHG4Hit, but we need to know which TrkrHitSet the TrkrHit is in
  // check if this association already exists
  // We need all hitsets with this key
  const auto hitsetrange = m_map.equal_range(hitsetkey);
  if( std::any_of(
    hitsetrange.first,
    hitsetrange.second,
    [&hitkey, &g4hitkey]( const MMap::value_type& pair ) { return pair.second.first == hitkey && pair.second.second == g4hitkey; } ) )
  { return; }

  // Does not exist, create it
  const auto assoc = std::make_pair(hitkey, g4hitkey);
  m_map.insert(hitsetrange.second, std::make_pair(hitsetkey, std::move(assoc)));
}

void
TrkrHitTruthAssocv1::removeAssoc(const TrkrDefs::hitsetkey hitsetkey, const TrkrDefs::hitkey hitkey)
{
  // remove all entries for this TrkrHit and its PHG4Hits, but we need to know which TrkrHitSet the TrkrHit is in
  // check if this association already exists
  // We need all hitsets with this key
  const auto hitsetrange = m_map.equal_range(hitsetkey);
  for( auto mapiter = hitsetrange.first; mapiter != hitsetrange.second; ++mapiter)
  {
    if(mapiter->second.first == hitkey)
    {
      m_map.erase (mapiter);
      return;
    }
  }

}

void
TrkrHitTruthAssocv1::getG4Hits(const TrkrDefs::hitsetkey hitsetkey, const unsigned int hidx, MMap &temp_map) const
{
  const auto hitsetrange = m_map.equal_range(hitsetkey);
  std::copy_if( hitsetrange.first, hitsetrange.second, std::inserter( temp_map, temp_map.end() ), 
    [hidx]( MMap::const_reference pair ) { return pair.second.first == hidx; } );
}
