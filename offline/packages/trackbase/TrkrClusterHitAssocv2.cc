/**
 * @file trackbase/TrkrClusterHitAssocv2.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief TrkrClusterHitAssoc implementation
 */

#include "TrkrClusterHitAssocv2.h"
#include "TrkrDefs.h"

#include <ostream>  // for operator<<, endl, basic_ostream, ostream, basic_o...

namespace
{
  TrkrClusterHitAssocv2::Map dummy_map;
}

//_________________________________________________________________________
void TrkrClusterHitAssocv2::Reset()
{
  for (unsigned int layer = 0; layer < max_layer; layer++)
  {
    for (unsigned int phi_segment = 0; phi_segment < max_phisegment; phi_segment++)
    {
      for (unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++)
      {
        m_map[layer][phi_segment][z_segment].clear();
      }
    }
  }
}

//_________________________________________________________________________
void TrkrClusterHitAssocv2::identify(std::ostream& os) const
{
  std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator iter;
  os << "-----TrkrClusterHitAssocv2-----" << std::endl;
  os << "Number of associations: " << size() << std::endl;
  for (unsigned int layer = 0; layer < max_layer; layer++)
  {
    for (unsigned int phi_segment = 0; phi_segment < max_phisegment; phi_segment++)
    {
      for (unsigned int z_segment = 0; z_segment < max_zsegment; z_segment++)
      {
        for (const auto& pair : m_map[layer][phi_segment][z_segment])
        {
          os << "clus key " << pair.first << std::dec
             << " layer " << TrkrDefs::getLayer(pair.first)
             << " hit key: " << pair.second << std::endl;
        }
      }
    }
  }
  os << "------------------------------" << std::endl;

  return;
}

//_________________________________________________________________________
void TrkrClusterHitAssocv2::addAssoc(TrkrDefs::cluskey ckey, unsigned int hidx)
{
  unsigned int layer = TrkrDefs::getLayer(ckey);
  unsigned int sector = TrkrDefs::getPhiElement(ckey);
  unsigned int side = TrkrDefs::getZElement(ckey);

  // bound check
  if (layer < max_layer && sector < max_phisegment && side < max_zsegment)
  {
    m_map[layer][sector][side].insert(std::make_pair(ckey, hidx));
  }
  else
  {
    std::cout
        << "TrkrClusterHitAssocv2::addAssoc - out of range access."
        << " layer: " << layer
        << " sector: " << sector
        << " side: " << side
        << std::endl;
  }
}

//_________________________________________________________________________
TrkrClusterHitAssocv2::Map* TrkrClusterHitAssocv2::getClusterMap(TrkrDefs::hitsetkey hitsetkey)
{
  unsigned int layer = TrkrDefs::getLayer(hitsetkey);
  unsigned int sector = TrkrDefs::getPhiElement(hitsetkey);
  unsigned int side = TrkrDefs::getZElement(hitsetkey);

  // bound check
  if (layer < max_layer && sector < max_phisegment && side < max_zsegment)
  {
    return &m_map[layer][sector][side];
  }
  else
  {
    std::cout
        << "TrkrClusterHitAssocv2::getClusterMap - out of range access."
        << " layer: " << layer
        << " sector: " << sector
        << " side: " << side
        << std::endl;
    return nullptr;
  }
}

//_________________________________________________________________________
TrkrClusterHitAssocv2::ConstRange TrkrClusterHitAssocv2::getHits(TrkrDefs::cluskey ckey)
{
  unsigned int layer = TrkrDefs::getLayer(ckey);
  unsigned int sector = TrkrDefs::getPhiElement(ckey);
  unsigned int side = TrkrDefs::getZElement(ckey);

  // bound check
  if (layer < max_layer && sector < max_phisegment && side < max_zsegment)
  {
    return std::make_pair(m_map[layer][sector][side].lower_bound(ckey), m_map[layer][sector][side].upper_bound(ckey));
  }
  else
  {
    std::cout
        << "TrkrClusterHitAssocv2::getHits - out of range access."
        << " layer: " << layer
        << " sector: " << sector
        << " side: " << side
        << std::endl;
    return std::make_pair(dummy_map.cbegin(), dummy_map.cend());
  }
}

//_________________________________________________________________________
unsigned int TrkrClusterHitAssocv2::size() const
{
  unsigned int size = 0;
  for (unsigned layer = 0; layer < max_layer; layer++)
  {
    for (unsigned phi_segment = 0; phi_segment < max_phisegment; phi_segment++)
    {
      for (unsigned z_segment = 0; z_segment < max_zsegment; z_segment++)
      {
        size += m_map[layer][phi_segment][z_segment].size();
      }
    }
  }

  return size;
}
