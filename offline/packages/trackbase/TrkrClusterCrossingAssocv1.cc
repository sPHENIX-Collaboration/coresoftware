/**
 * @file trackbase/TrkrClusterCrossingAssocv1.cc
 * @author Tony Frawley
 * @date March 2022
 * @brief TrkrClusterCrossingAssoc implementation
 */

#include "TrkrClusterCrossingAssocv1.h"
#include "TrkrDefs.h"

#include <ostream>  // for operator<<, endl, basic_ostream, ostream, basic_o...

namespace
{
  TrkrClusterCrossingAssocv1::Map dummy_map;
}

//_________________________________________________________________________
void TrkrClusterCrossingAssocv1::Reset()
{ 
 
 
 // delete all entries
    Map empty_map;
    empty_map.swap(m_map);
   empty_map.clear();

   m_map.clear();
}

//_________________________________________________________________________
void TrkrClusterCrossingAssocv1::identify(std::ostream &os) const
{
  std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator iter;
  os << "-----TrkrClusterCrossingAssocv1-----" << std::endl;
  os << "Number of associations: " << size() << std::endl;
  for( const auto& map_pair:m_map )
    {
      os << "clus key " << map_pair.first << std::dec
        << " layer " << (unsigned int) TrkrDefs::getLayer(map_pair.first)
        << " crossing number: " << map_pair.second << std::endl;
    }
  os << "------------------------------" << std::endl;

  return;
  
}

//_________________________________________________________________________
void TrkrClusterCrossingAssocv1::addAssoc(TrkrDefs::cluskey ckey, short int hidx)
{
  // insert association
  m_map.insert(std::make_pair(ckey, hidx));
}

//_________________________________________________________________________
TrkrClusterCrossingAssocv1::ConstRange TrkrClusterCrossingAssocv1::getCrossings(TrkrDefs::cluskey ckey) const
{
  const auto range = m_map.equal_range(ckey);
  return range;    
}

//_________________________________________________________________________
unsigned int TrkrClusterCrossingAssocv1::size(void) const
{
  unsigned int size = 0;
  size = m_map.size();

  return size;
}

TrkrClusterCrossingAssocv1::ConstRange TrkrClusterCrossingAssocv1::getAll() const
{
  return std::make_pair(m_map.begin(), m_map.end());
}
