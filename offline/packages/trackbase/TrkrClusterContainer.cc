/**
 * @file trackbase/TrkrClusterContainer.cc
 * @author D. McGlinchey, Hugo Pereira Da Costa
 * @date June 2018
 */
#include "TrkrClusterContainer.h"

namespace
{
  TrkrClusterContainer::Map dummy_map;
}

//__________________________________________________________
TrkrClusterContainer::ConstRange TrkrClusterContainer::getClusters() const
{
  return std::make_pair(dummy_map.cbegin(), dummy_map.cend());
}

//__________________________________________________________
TrkrClusterContainer::ConstRange TrkrClusterContainer::getClusters(TrkrDefs::hitsetkey /*hitsetkey*/)
{
  return std::make_pair(dummy_map.cbegin(), dummy_map.cend());
}
