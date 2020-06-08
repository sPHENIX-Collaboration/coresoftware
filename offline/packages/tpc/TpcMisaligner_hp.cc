#include "TpcMisaligner_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <math.h>

//_____________________________________________________________________
template<class T> constexpr T square( const T& x ) { return x*x; }

//_____________________________________________________________________
TpcMisaligner_hp::TpcMisaligner_hp( const std::string& name ):
  SubsysReco( name)
  {}

//_____________________________________________________________________
int TpcMisaligner_hp::Init(PHCompositeNode*)
{ return Fun4AllReturnCodes::EVENT_OK; }

//_____________________________________________________________________
int TpcMisaligner_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  transform_clusters();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcMisaligner_hp::load_nodes( PHCompositeNode* topNode )
{
  // cluster map
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcMisaligner_hp::transform_clusters()
{
  if( !_cluster_map ) return;

  auto range = _cluster_map->getClusters();
  for( auto clusterIter = range.first; clusterIter != range.second; ++clusterIter )
  {
    // check if cluster belongs to TPC
    const auto& key = clusterIter->first;
    const auto trkrid = TrkrDefs::getTrkrId(key);
    if( trkrid == TrkrDefs::tpcId )
    { transform_cluster( clusterIter->second ); }
  }

  return;
}

//_____________________________________________________________________
void TpcMisaligner_hp::transform_cluster( TrkrCluster* cluster )
{

  // define the r transformation (cm)
  static constexpr float rmin_tpc = 20;
  static constexpr float rmax_tpc = 78;
  static constexpr float rlength_tpc = rmax_tpc - rmin_tpc;

  const float r = std::sqrt( square( cluster->getX() ) + square( cluster->getY() ) );
  const float phi = std::atan2( cluster->getY(), cluster->getX() );
  const float z = cluster->getZ();

  #if true
  // rphi distortion
  static constexpr float deltarphi_max = 1.5;
  const float deltarphi = deltarphi_max*std::cos(2.*M_PI*(r-rmin_tpc)/rlength_tpc);
  const float deltaphi = deltarphi/r;
  #else
  const float deltaphi = 0;
  #endif

  #if false
  // z distortion
  static constexpr float deltaz_max = 1.5;
  const float deltaz = deltaz_max*std::cos(2.*M_PI*(r-rmin_tpc)/rlength_tpc);
  #else
  const float deltaz = 0;
  #endif

  #if false
  // r distortion
  static constexpr float deltar_max = 1.5;
  const float deltar = deltar_max*std::cos(2.*M_PI*(r-rmin_tpc)/rlength_tpc);
  #else
  const float deltar = 0;
  #endif

  // modify cluster x,y positions
  const float newr = r+deltar;
  const float newphi = phi + deltaphi;
  const float cosnewphi = std::cos( newphi );
  const float sinnewphi = std::sin( newphi );
  const float newx = newr*cosnewphi;
  const float newy = newr*sinnewphi;
  const float newz = z + deltaz;

  cluster->setX( newx );
  cluster->setY( newy );
  cluster->setZ( newz );

}
