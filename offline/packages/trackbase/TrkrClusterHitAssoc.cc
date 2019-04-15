/**                                                                                                                                                                                                       
 * @file trackbase/TrkrClusterHitAssoc.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief TrkrClusterHitAssoc implementation
 */

#include "TrkrClusterHitAssoc.h"

#include <cstdlib>

TrkrClusterHitAssoc::TrkrClusterHitAssoc() 
  : m_map()
{
}


TrkrClusterHitAssoc::~TrkrClusterHitAssoc()
{
}

void 
TrkrClusterHitAssoc::Reset()
{
  m_map.clear();
}

void 
TrkrClusterHitAssoc::identify(std::ostream &os) const
{
  os << "-----TrkrClusterHitAssoc-----" << std::endl;
  os << "Number of associations: " << m_map.size() << std::endl;

  for ( auto& entry : m_map )
  {
    // os << "   cluster key: 0x" << std::hex << entry.first << std::dec
    os << "   cluster key: "  << entry.first << std::dec
       << " hit key: " << entry.second << std::endl;
  }

  os << "------------------------------" << std::endl;

  return;

}

void 
TrkrClusterHitAssoc::addAssoc(TrkrDefs::cluskey ckey, unsigned int hidx)
{
  m_map.insert(std::make_pair(ckey, hidx));
}

TrkrClusterHitAssoc::ConstRange 
TrkrClusterHitAssoc::getHits(TrkrDefs::cluskey ckey)
{
  ConstRange retpair;
  retpair.first = m_map.lower_bound(ckey);
  retpair.second = m_map.upper_bound(ckey);
  return retpair;
}
