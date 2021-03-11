#include "TpcSpaceChargeCorrection.h"

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

  const std::array<const std::string,2> extension = {{ "_negz", "_posz" }};
  for( int i =0; i < 2; ++i )
  {
    m_hDPint[i] = dynamic_cast<TH3*>(m_distortion_tfile->Get(Form("hIntDistortionP%s", extension[i].c_str()))); assert( m_hDPint[i] );
    m_hDRint[i] = dynamic_cast<TH3*>(m_distortion_tfile->Get(Form("hIntDistortionR%s", extension[i].c_str()))); assert( m_hDRint[i] );
    m_hDZint[i] = dynamic_cast<TH3*>(m_distortion_tfile->Get(Form("hIntDistortionZ%s", extension[i].c_str()))); assert( m_hDZint[i] );
  }

  // coordinates
  std::cout << "TpcSpaceChargeCorrection::InitRun - coordinates: " << std::bitset<3>(m_coordinates) << std::endl;

  // dump axis limits
  if( Verbosity() )
  {
    for( int i =0; i < 2; ++i )
    {
      std::cout << "TpcSpaceChargeCorrection::InitRun - histogram: " << m_hDPint[i]->GetName() << std::endl;
      for(const auto& axis:{ m_hDPint[i]->GetXaxis(), m_hDPint[i]->GetYaxis(), m_hDPint[i]->GetZaxis() })
      {
        std::cout
          << "TpcSpaceChargeCorrection::InitRun -"
          << " axis: " << axis->GetTitle()
          << " bins: " << axis->GetNbins()
          << " limits: " << axis->GetXmin() << " " << axis->GetXmax()
          << std::endl;
      }
    }
  }

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
  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void TpcSpaceChargeCorrection::transform_clusters()
{
  if( !m_cluster_map ) return;

  auto range = m_cluster_map->getClusters();
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
void TpcSpaceChargeCorrection::transform_cluster( TrkrCluster* cluster )
{
  // get cluster radius, phi and z
  const auto r = std::sqrt( square( cluster->getX() ) + square( cluster->getY() ) );
  auto phi = std::atan2( cluster->getY(), cluster->getX() );
  if( phi < 0 ) phi += 2*M_PI;

  const auto z = cluster->getZ();
  const int index = z > 0 ? 1:0;

  // apply corrections
  const auto phi_new = (m_coordinates & COORD_PHI) ? phi - m_hDPint[index]->Interpolate(phi,r,z)/r : phi;
  const auto r_new = (m_coordinates & COORD_R) ? r - m_hDRint[index]->Interpolate(phi,r,z) : r;
  const auto z_new = (m_coordinates & COORD_Z) ? z - m_hDZint[index]->Interpolate(phi,r,z) : z;

  // update cluster
  const auto x_new = r_new*std::cos( phi_new );
  const auto y_new = r_new*std::sin( phi_new );

  cluster->setX( x_new );
  cluster->setY( y_new );
  cluster->setZ( z_new );
}
