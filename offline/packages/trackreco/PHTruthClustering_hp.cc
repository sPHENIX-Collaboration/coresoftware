#include "PHTruthClustering_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrHitTruthAssoc.h>

//_____________________________________________________________________
namespace
{

  /// square
  template<class T> T square( T x ) { return x*x; }

  /// radius
  template<class T> T get_r( T x, T y ) { return std::sqrt( square(x) + square(y) ); }

  /// radius
  float get_r( PHG4Hit* hit, int i )
  {  return get_r( hit->get_x(i), hit->get_y(i) ); }

  /// calculate the average of member function called on all members in collection
  template< float (PHG4Hit::*accessor)(int) const>
  float interpolate( std::set<PHG4Hit*> hits, float rextrap )
  {
    // calculate all terms needed for the interpolation
    // need to use double everywhere here due to numerical divergences
    double sw = 0;
    double swr = 0;
    double swr2 = 0;
    double swx = 0;
    double swrx = 0;

    bool valid( false );
    for( const auto& hit:hits )
    {

      const double x0 = (hit->*accessor)(0);
      const double x1 = (hit->*accessor)(1);
      if( std::isnan( x0 ) || std::isnan( x1 ) ) continue;

      const double w = hit->get_edep();
      if( w < 0 ) continue;

      valid = true;
      const double r0 = get_r( hit, 0 );
      const double r1 = get_r( hit, 1 );

      sw += w*2;
      swr += w*(r0 + r1);
      swr2 += w*(square(r0) + square(r1));
      swx += w*(x0 + x1);
      swrx += w*(r0*x0 + r1*x1);
    }

    if( !valid ) return NAN;

    const auto alpha = (sw*swrx - swr*swx);
    const auto beta = (swr2*swx - swr*swrx);
    const auto denom = (sw*swr2 - square(swr));

    return ( alpha*rextrap + beta )/denom;
  }

}

//_____________________________________________________________________
PHTruthClustering_hp::PHTruthClustering_hp( const std::string& name ):
  SubsysReco( name)
{}

//_____________________________________________________________________
int PHTruthClustering_hp::Init(PHCompositeNode* topNode )
{

  // find DST node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "PHTruthClustering_hp::Init - DST Node missing" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int PHTruthClustering_hp::InitRun(PHCompositeNode* )
{
  std::cout << "PHTruthClustering_hp::InitRun." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int PHTruthClustering_hp::process_event(PHCompositeNode* topNode)
{

  // load nodes
  auto res = load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  replace_clusters();

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int PHTruthClustering_hp::End(PHCompositeNode* )
{
  std::cout << "PHTruthClustering_hp::End." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int PHTruthClustering_hp::load_nodes( PHCompositeNode* topNode )
{

  // cluster map
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  // cluster hit association map
  _cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");

  // cluster hit association map
  _hit_truth_map = findNode::getClass<TrkrHitTruthAssoc>(topNode,"TRKR_HITTRUTHASSOC");

  // g4hits
  _g4hits_tpc = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  _g4hits_intt = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_INTT");
  _g4hits_mvtx = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MVTX");
  _g4hits_micromegas = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MICROMEGAS");

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
void PHTruthClustering_hp::replace_clusters()
{

  if( !_cluster_map ) return;

  auto range = _cluster_map->getClusters();
  for( auto clusterIter = range.first; clusterIter != range.second; ++clusterIter )
  {
    const auto& key = clusterIter->first;
    const auto& cluster = clusterIter->second;
    const auto g4hits = find_g4hits( key );
    if( g4hits.empty() ) continue;

    // get averaged coordinates and assign to cluster
    const auto rextrap = get_r( cluster->getX(), cluster->getY() );
    cluster->setX( interpolate<&PHG4Hit::get_x>( g4hits, rextrap ) );
    cluster->setY( interpolate<&PHG4Hit::get_y>( g4hits, rextrap ) );
    cluster->setZ( interpolate<&PHG4Hit::get_z>( g4hits, rextrap ) );
  }

}

//_____________________________________________________________________
PHTruthClustering_hp::G4HitSet PHTruthClustering_hp::find_g4hits( TrkrDefs::cluskey cluster_key ) const
{

  // check maps
  if( !( _cluster_hit_map && _hit_truth_map ) ) return G4HitSet();

  // find hitset associated to cluster
  G4HitSet out;
  const auto hitset_key = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);

  // loop over hits associated to clusters
  const auto range = _cluster_hit_map->getHits(cluster_key);
  for( auto iter = range.first; iter != range.second; ++iter )
  {

    // hit key
    const auto& hit_key = iter->second;

    // store hits to g4hit associations
    TrkrHitTruthAssoc::MMap g4hit_map;
    _hit_truth_map->getG4Hits( hitset_key, hit_key, g4hit_map );

    // find corresponding g4 hist
    for( auto truth_iter = g4hit_map.begin(); truth_iter != g4hit_map.end(); ++truth_iter )
    {

      const auto g4hit_key = truth_iter->second.second;
      PHG4Hit* g4hit = nullptr;

      switch( TrkrDefs::getTrkrId( hitset_key ) )
      {
        case TrkrDefs::tpcId:
        if( _g4hits_tpc ) g4hit = _g4hits_tpc->findHit( g4hit_key );
        break;

        case TrkrDefs::inttId:
        if( _g4hits_intt ) g4hit = _g4hits_intt->findHit( g4hit_key );
        break;

        case TrkrDefs::micromegasId:
        if( _g4hits_micromegas ) g4hit = _g4hits_micromegas->findHit( g4hit_key );
        break;

        case TrkrDefs::mvtxId:
        if( _g4hits_mvtx ) g4hit = _g4hits_mvtx->findHit( g4hit_key );
        break;

        default: break;
      }

      if( g4hit ) out.insert( g4hit );
      else std::cout << "PHTruthClustering_hp::find_g4hits - g4hit not found " << g4hit_key << std::endl;

    }
  }

  // insert in map and return
  return out;

}
