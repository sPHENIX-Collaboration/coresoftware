/**
 * @file trackbase/TrkrClusterContainer.cc
 * @author D. McGlinchey Hugo PEREIRA DA COSTA
 * @date June 2018
 */
#include "TrkrClusterContainer.h"

namespace
{
  TrkrClusterContainer::Map dummy_map;
}

//__________________________________________________________
ConstIterator addCluster(TrkrCluster*)
{ return dummy_map.cbegin(); }

//__________________________________________________________
ConstIterator addClusterSpecifyKey(const TrkrDefs::cluskey, TrkrCluster* )
{ return dummy_map.cbegin(); }

//__________________________________________________________
Iterator findOrAddCluster(TrkrDefs::cluskey)
{ return dummy_map.begin(); }

//__________________________________________________________
ConstRange getClusters(TrkrDefs::hitsetkey hitsetkey) const
{ return std::make_pair( dummy_map.cbegin(), dummy_map.cend() ); }
