#include "TpcSpaceChargeCorrection.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
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
TpcSpaceChargeCorrection::TpcSpaceChargeCorrection( const std::string& name ):
  SubsysReco( name)
  {}

//_____________________________________________________________________
int TpcSpaceChargeCorrection::InitRun(PHCompositeNode*)
{

  std::cout << "TpcSpaceChargeCorrection::InitRun - reading distortions from " << m_distortion_filename << std::endl;
  m_distortion_tfile = TFile::Open( m_distortion_filename.c_str());
  if( !m_distortion_tfile )
  {
    std::cout << "TpcSpaceChargeCorrection::InitRun - cannot open " << m_distortion_filename << std::endl;
    exit(1);
  }


  // Open TH3F files only once that contain distortions due to space charge
  m_hDPint= dynamic_cast<TH3*>(m_distortion_tfile->Get("hIntDistortionP")); assert( m_hDPint );
  m_hDRint= dynamic_cast<TH3*>(m_distortion_tfile->Get("hIntDistortionR")); assert( m_hDRint );
  m_hDZint= dynamic_cast<TH3*>(m_distortion_tfile->Get("hIntDistortionZ")); assert( m_hDZint );

  // coordinates
  std::cout << "TpcSpaceChargeCorrection::InitRun - coordinates: " << std::bitset<3>(m_coordinates) << std::endl;

  m_fullzrange = (m_hDPint->GetZaxis()->GetXmin() == -m_hDPint->GetZaxis()->GetXmax());
  std::cout << "TpcSpaceChargeCorrection::InitRun - m_fullzrange = " << std::boolalpha << m_fullzrange << std::endl;

  // dump axis limits
  for(const auto& axis:{ m_hDPint->GetXaxis(), m_hDPint->GetYaxis(), m_hDPint->GetZaxis() })
  { std::cout << "TpcSpaceChargeCorrection::InitRun - axis: " << axis->GetTitle() << " bins: " << axis->GetNbins() << " limits: " << axis->GetXmin() << " " << axis->GetXmax() << std::endl; }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeCorrection::process_event(PHCompositeNode* topNode)
{
  // load nodes
  auto res =  load_nodes(topNode);
  if( res != Fun4AllReturnCodes::EVENT_OK ) return res;

  transform_clusters();
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcSpaceChargeCorrection::load_nodes( PHCompositeNode* topNode )
{
  // cluster map
  m_hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcSpaceChargeCorrection::transform_clusters()
{
  if( !m_cluster_map ) return;
  if( !m_hitsets ) return;
  
  auto hitsetrange = m_hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (auto hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr){
    auto range = m_cluster_map->getClusters(hitsetitr->first);
    for( auto clusterIter = range.first; clusterIter != range.second; ++clusterIter )
      {
	transform_cluster( clusterIter->second );
      }
  }

  return;
}

//_____________________________________________________________________
void TpcSpaceChargeCorrection::transform_cluster( TrkrCluster* cluster )
{
  // get cluster radius, phi and z
  const auto r = std::sqrt( square( cluster->getX() ) + square( cluster->getY() ) );
  auto phi = std::atan2( cluster->getY(), cluster->getX() );
  if( phi < 0 ) phi += 2*M_PI;

  const auto z = cluster->getZ();
  const auto zmap = m_fullzrange ? z:std::abs(z);

  // apply corrections
  const auto phi_new = (m_coordinates & COORD_PHI) ? phi - m_hDPint->Interpolate(phi,r,zmap)/r : phi;
  const auto r_new = (m_coordinates & COORD_R) ? r - m_hDRint->Interpolate(phi,r,zmap) : r;
  const auto z_new = (m_coordinates & COORD_Z) ? z - m_hDZint->Interpolate(phi,r,zmap) : z;

  // update cluster
  const auto x_new = r_new*std::cos( phi_new );
  const auto y_new = r_new*std::sin( phi_new );

  cluster->setX( x_new );
  cluster->setY( y_new );
  cluster->setZ( z_new );
}
