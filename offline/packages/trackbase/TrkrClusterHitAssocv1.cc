/**                                                                                                                                                                                                       
 * @file trackbase/TrkrClusterHitAssoc.cc
 * @author D. McGlinchey
 * @date June 2018
 * @brief TrkrClusterHitAssoc implementation
 */

#include "TrkrClusterHitAssocv1.h"

#include <ostream>  // for operator<<, endl, basic_ostream, ostream, basic_o...

TrkrClusterHitAssocv1::TrkrClusterHitAssocv1() 
  : m_map()
{
}


TrkrClusterHitAssocv1::~TrkrClusterHitAssocv1()
{
}

void 
TrkrClusterHitAssocv1::Reset()
{
  m_map.clear();
}

void 
TrkrClusterHitAssocv1::identify(std::ostream &os) const
{
  os << "-----TrkrClusterHitAssocv1-----" << std::endl;
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
TrkrClusterHitAssocv1::addAssoc(TrkrDefs::cluskey ckey, unsigned int hidx)
{
  m_map.insert(std::make_pair(ckey, hidx));
}

std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>
TrkrClusterHitAssocv1::getHits(TrkrDefs::cluskey ckey)
{
  std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>  retpair;
  retpair.first = m_map.lower_bound(ckey);
  retpair.second = m_map.upper_bound(ckey);
  return retpair;
}
