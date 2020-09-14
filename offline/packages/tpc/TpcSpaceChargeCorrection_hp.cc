#include "TpcSpaceChargeCorrection_hp.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <TFile.h>
#include <TH3.h>

#include <bitset>
#include <math.h>

namespace
{
  //_____________________________________________________________________
  template<class T> constexpr T square( const T& x ) { return x*x; }
}

//_____________________________________________________________________
TpcSpaceChargeCorrection_hp::TpcSpaceChargeCorrection_hp( const std::string& name ):
  SubsysReco( name)
  {}

//_____________________________________________________________________
int TpcSpaceChargeCorrection_hp::InitRun(PHCompositeNode*)
{

  std::cout << "TpcSpaceChargeCorrection_hp::InitRun - reading distortions from " << m_distortion_filename << std::endl;
  m_distortion_tfile = TFile::Open( m_distortion_filename.c_str());
  if( !m_distortion_tfile )
  {
    std::cout << "TpcSpaceChargeCorrection_hp::InitRun - cannot open " << m_distortion_filename << std::endl;
    exit(1);
  }


  // Open TH3F files only once that contain distortions due to space charge
  hDPint= dynamic_cast<TH3*>(m_distortion_tfile->Get("hIntDistortionP")); assert( hDPint );
  hDRint= dynamic_cast<TH3*>(m_distortion_tfile->Get("hIntDistortionR")); assert( hDRint );
  hDZint= dynamic_cast<TH3*>(m_distortion_tfile->Get("hIntDistortionZ")); assert( hDZint );

  // coordinates
  std::cout << "TpcSpaceChargeCorrection_hp::InitRun - coordinates: " << std::bitset<3>(m_coordinates) << std::endl;

  m_fullzrange = (hDPint->GetZaxis()->GetXmin() == -hDPint->GetZaxis()->GetXmax());
  std::cout << "TpcSpaceChargeCorrection_hp::InitRun - m_fullzrange = " << m_fullzrange << std::endl;

  // dump axis limits
  for(const auto& axis:{ hDPint->GetXaxis(), hDPint->GetYaxis(), hDPint->GetZaxis() })
  { std::cout << "TpcSpaceChargeCorrection_hp::InitRun - axis: " << axis->GetTitle() << " bins: " << axis->GetNbins() << " limits: " << axis->GetXmin() << " " << axis->GetXmax() << std::endl; }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeCorrection_hp::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  transform_clusters();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeCorrection_hp::load_nodes( PHCompositeNode* topNode )
{
  // cluster map
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcSpaceChargeCorrection_hp::transform_clusters()
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
void TpcSpaceChargeCorrection_hp::transform_cluster( TrkrCluster* cluster )
{
  // get cluster radius, phi and z
  const auto r = std::sqrt( square( cluster->getX() ) + square( cluster->getY() ) );
  auto phi = std::atan2( cluster->getY(), cluster->getX() );
  if( phi < 0 ) phi += 2*M_PI;

  const auto z = cluster->getZ();
  const auto zmap = m_fullzrange ? z:std::abs(z);

  // apply corrections
  const auto phi_new = (m_coordinates & COORD_PHI) ? phi - hDPint->Interpolate(phi,r,zmap)/r : phi;

  // dont move r for the moment
  const auto r_new = (m_coordinates & COORD_R) ? r - hDRint->Interpolate(phi,r,zmap) : r;

  // dont move z for the moment
  const auto z_new = (m_coordinates & COORD_R) ? z - hDZint->Interpolate(phi,r,zmap) : z;

  // update cluster
  const auto x_new = r_new*std::cos( phi_new );
  const auto y_new = r_new*std::sin( phi_new );

  cluster->setX( x_new );
  cluster->setY( y_new );
  cluster->setZ( z_new );
}
