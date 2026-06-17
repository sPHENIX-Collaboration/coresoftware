#include "PHTrackTrackSeedSynchronization.h"

/// Tracking includes
#include <trackbase/MvtxDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrDefs.h>  // for cluskey, getTrkrId, tpcId

#include <trackbase_historic/SvtxTrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/sphenix_constants.h>

#include <TF1.h>
#include <TFile.h>
#include <TNtuple.h>

#include <climits>   // for UINT_MAX
#include <cmath>     // for fabs, sqrt
#include <iostream>  // for operator<<, basic_ostream
#include <memory>
#include <set>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

using namespace std;

//____________________________________________________________________________..
PHTrackTrackSeedSynchronization::PHTrackTrackSeedSynchronization(const std::string &name)
  : SubsysReco(name)
{}

//____________________________________________________________________________..
int PHTrackTrackSeedSynchronization::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  return ret;
}

//____________________________________________________________________________..
int PHTrackTrackSeedSynchronization::process_event(PHCompositeNode * /*unused*/)
{
  // loop over tracks and synchronize
  for( auto &&[key,track]:*_svtx_track_map )
  { synchronize_track(track); }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________
bool PHTrackTrackSeedSynchronization::synchronize_track( SvtxTrack* track ) const
{
  {
    // silicon seed
    auto* seed = track->get_silicon_seed();
    if( seed )
    {
      auto index = find_seed_id( _si_seed_map, seed );
      if( index < _si_seed_map->size() )
      { track->set_silicon_seed( _si_seed_map->get(index)); }
    }
  }

  {
    // tpc seed
    auto* seed = track->get_tpc_seed();
    if( seed )
    {
      auto index = find_seed_id( _tpc_seed_map, seed );
      if( index < _tpc_seed_map->size() )
      { track->set_tpc_seed( _tpc_seed_map->get(index)); }
    }
  }

  return true;
}

//__________________________________________________________________________________
int PHTrackTrackSeedSynchronization::End(PHCompositeNode * /*unused*/)
{ return Fun4AllReturnCodes::EVENT_OK; }

//__________________________________________________________________________________
int PHTrackTrackSeedSynchronization::GetNodes(PHCompositeNode *topNode)
{

  // tracks
  _svtx_track_map = findNode::getClass<SvtxTrackMap>(topNode, _svtx_track_map_name);
  if (!_svtx_track_map)
  {
    cerr << "PHTrackTrackSeedSynchronization::GetNodes - " << _svtx_track_map_name.c_str() << " not found on the node tree." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // silicon seeds
  _si_seed_map = findNode::getClass<TrackSeedContainer>(topNode, _si_seed_map_name);
  if (!_si_seed_map)
  {
    cerr << "PHTrackTrackSeedSynchronization::GetNodes - " << _svtx_track_map_name.c_str() << " not found on the node tree." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc seeds
  _tpc_seed_map = findNode::getClass<TrackSeedContainer>(topNode, _tpc_seed_map_name);
  if (!_tpc_seed_map)
  {
    cerr << "PHTrackTrackSeedSynchronization::GetNodes - " << _svtx_track_map_name.c_str() << " not found on the node tree." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


//__________________________________________________________________________________
size_t PHTrackTrackSeedSynchronization::find_seed_id( TrackSeedContainer* container, TrackSeed* source ) const
{
  // perform quick search
  const size_t index = container->find( source );
  if( index < container->size() ) return index;

  // perform deep search based on cluster keys
  if( Verbosity() )
  { std::cout << "PHTrackTrackSeedSynchronization::find_seed_id - performing deep search for seed " << source << " in container " << container << std::endl; }

  using cluster_keyset_t=std::set<TrkrDefs::cluskey>;

  // get cluster key set from seed
  auto get_cluster_keyset = [this]( TrackSeed* seed )
  {
    cluster_keyset_t ckeys;
    if( m_ignore_micromegas )
    {
      std::copy_if( seed->begin_cluster_keys(), seed->end_cluster_keys(), std::inserter(ckeys, ckeys.end() ),
        []( const TrkrDefs::cluskey& ckey ) { return TrkrDefs::getTrkrId(ckey) != TrkrDefs::micromegasId; } );
    } else {
      std::copy( seed->begin_cluster_keys(), seed->end_cluster_keys(), std::inserter(ckeys, ckeys.end() ) );
    }
    return ckeys;

  };

  const auto source_ckeys = get_cluster_keyset( source );
  for( size_t i = 0; i < container->size(); ++i )
  {
    auto* seed = container->get(i);
    if( !seed ) continue;

    const auto ckeys = get_cluster_keyset( seed );
    if( ckeys == source_ckeys )
    { return i; }
  }

  // error
  std::cout << "PHTrackTrackSeedSynchronization::find_seed_id - could not find seed " << source << " in container " << container << std::endl;
  return container->size();
}
