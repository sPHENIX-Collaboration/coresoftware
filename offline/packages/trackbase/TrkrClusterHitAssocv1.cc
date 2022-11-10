/**
 * @file trackbase/TrkrClusterHitAssoc.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief TrkrClusterHitAssoc implementation
 */

#include "TrkrClusterHitAssocv1.h"

#include <ostream>  // for operator<<, endl, basic_ostream, ostream, basic_o...

void TrkrClusterHitAssocv1::Reset()
{
  m_map.clear();
}

void TrkrClusterHitAssocv1::identify(std::ostream& os) const
{
  os << "-----TrkrClusterHitAssocv1-----" << std::endl;
  os << "Number of associations: " << m_map.size() << std::endl;

  for (auto& entry : m_map)
  {
    // os << "   cluster key: 0x" << std::hex << entry.first << std::dec
    os << "   cluster key: " << entry.first << std::dec
       << " hit key: " << entry.second << std::endl;
  }

  os << "------------------------------" << std::endl;

  return;
}

void TrkrClusterHitAssocv1::addAssoc(TrkrDefs::cluskey ckey, unsigned int hidx)
{
  m_map.insert(std::make_pair(ckey, hidx));
}

TrkrClusterHitAssocv1::ConstRange
TrkrClusterHitAssocv1::getHits(TrkrDefs::cluskey ckey)
{
  return std::make_pair(m_map.lower_bound(ckey), m_map.upper_bound(ckey));
}
