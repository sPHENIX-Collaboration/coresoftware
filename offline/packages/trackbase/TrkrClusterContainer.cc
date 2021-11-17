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
TrkrClusterContainer::ConstIterator TrkrClusterContainer::addCluster(TrkrCluster*)
{ return dummy_map.cbegin(); }

//__________________________________________________________
TrkrClusterContainer::ConstIterator TrkrClusterContainer::addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster* )
{ return dummy_map.cbegin(); }

//__________________________________________________________
TrkrClusterContainer::Iterator TrkrClusterContainer::findOrAddCluster(TrkrDefs::cluskey)
{ return dummy_map.begin(); }

//__________________________________________________________
TrkrClusterContainer::ConstRange TrkrClusterContainer::getClusters() const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }

//__________________________________________________________
TrkrClusterContainer::ConstRange TrkrClusterContainer::getClusters(TrkrDefs::hitsetkey /*hitsetkey*/) const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }
