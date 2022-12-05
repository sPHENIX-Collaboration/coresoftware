#include "ParticleFlowReco.h"

#include "ParticleFlowElementContainer.h"
#include "ParticleFlowElementv1.h"

#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawClusterUtility.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackState.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>


#include <TLorentzVector.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

#include <cmath>
#include <iostream>

// examine second value of std::pair, sort by smallest
bool sort_by_pair_second_lowest( const std::pair<int,float> &a,  const std::pair<int,float> &b) 
{ 
  return (a.second < b.second); 
} 

float ParticleFlowReco::calculate_dR( float eta1, float eta2, float phi1, float phi2 ) {

  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  while ( dphi > M_PI ) dphi -= 2 * M_PI;
  while ( dphi < -M_PI ) dphi += 2 * M_PI;
  return sqrt( pow( deta, 2 ) + pow( dphi ,2 ) );

}

std::pair<float, float> ParticleFlowReco::get_expected_signature( int trk ) {
  
  float response = ( 0.553437 + 0.0572246 * log( _pflow_TRK_p[ trk ] ) ) * _pflow_TRK_p[ trk ];
  float resolution = sqrt( pow( 0.119123 , 2 ) + pow( 0.312361 , 2 ) / _pflow_TRK_p[ trk ] ) * _pflow_TRK_p[ trk ] ;

  std::pair<float, float> expected_signature( response , resolution );

  return expected_signature;

}

//____________________________________________________________________________..
ParticleFlowReco::ParticleFlowReco(const std::string &name):
  SubsysReco(name),
  _energy_match_Nsigma( 1.5 )
{
}

//____________________________________________________________________________..
ParticleFlowReco::~ParticleFlowReco()
{
}

//____________________________________________________________________________..
int ParticleFlowReco::Init(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ParticleFlowReco::InitRun(PHCompositeNode *topNode)
{
  return CreateNode(topNode);

}

//____________________________________________________________________________..
int ParticleFlowReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
  std::cout << "ParticleFlowReco::process_event with Nsigma = " << _energy_match_Nsigma << std::endl;
  }
  // get handle to pflow node
  ParticleFlowElementContainer *pflowContainer = findNode::getClass<ParticleFlowElementContainer>(topNode, "ParticleFlowElements");
  if (!pflowContainer) {
    std::cout << " ERROR -- can't find ParticleFlowElements node after it should have been created" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // used for indexing PFlow elements in container
  int global_pflow_index = 0;
  
  // read in towers 
  RawTowerContainer *towersEM = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  RawTowerContainer *towersIH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");

  if ( !towersEM || !towersIH || !towersOH ) {
    std::cout << "ParticleFlowReco::process_event : FATAL ERROR, cannot find tower containers" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // read in tower geometries
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  if ( !geomEM || !geomIH || !geomOH ) {
    std::cout << "ParticleFlowReco::process_event : FATAL ERROR, cannot find tower geometry containers" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // read in clusters 
  RawClusterContainer *clustersEM = findNode::getClass<RawClusterContainer>(topNode, "TOPOCLUSTER_EMCAL");
  RawClusterContainer *clustersHAD = findNode::getClass<RawClusterContainer>(topNode, "TOPOCLUSTER_HCAL");

  if ( !clustersEM ) {
    std::cout << "ParticleFlowReco::process_event : FATAL ERROR, cannot find cluster container TOPOCLUSTER_EMCAL" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if ( !clustersHAD ) {
    std::cout << "ParticleFlowReco::process_event : FATAL ERROR, cannot find cluster container TOPOCLUSTER_HCAL" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // reset internal particle-flow representation
  _pflow_TRK_p.clear();
  _pflow_TRK_eta.clear();
  _pflow_TRK_phi.clear();
  _pflow_TRK_match_EM.clear();
  _pflow_TRK_match_HAD.clear();
  _pflow_TRK_addtl_match_EM.clear();
  _pflow_TRK_trk.clear();
  _pflow_TRK_EMproj_phi.clear();
  _pflow_TRK_EMproj_eta.clear();
  _pflow_TRK_HADproj_eta.clear();
  _pflow_TRK_HADproj_phi.clear();

  _pflow_EM_E.clear();
  _pflow_EM_eta.clear();
  _pflow_EM_phi.clear();
  _pflow_EM_tower_eta.clear();
  _pflow_EM_tower_phi.clear();
  _pflow_EM_match_HAD.clear();
  _pflow_EM_match_TRK.clear();
  _pflow_EM_cluster.clear();

  _pflow_HAD_E.clear();
  _pflow_HAD_eta.clear();
  _pflow_HAD_phi.clear();
  _pflow_HAD_tower_eta.clear();
  _pflow_HAD_tower_phi.clear();
  _pflow_HAD_match_EM.clear();
  _pflow_HAD_match_TRK.clear();
  _pflow_HAD_cluster.clear();

  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  GlobalVertex* vertex = nullptr;

  if(vertexmap)
    {
      if (!vertexmap->empty())
        {
	  vertex = (vertexmap->begin()->second);
	}
    }

  if ( Verbosity() > 2 ) 
    std::cout << "ParticleFlowReco::process_event : initial population of TRK, EM, HAD objects " << std::endl;

  // read in tracks with > 0.5 GeV
  {
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, _track_map_name);
  
    float cemcradius = geomEM->get_radius();
    float ohcalradius = geomOH->get_radius();

    for(auto iter = trackmap->begin(); iter != trackmap->end(); ++iter)
      {
	SvtxTrack *track = iter->second;

	if(track->get_pt() < 0.5)
	  { continue; }

	if(fabs(track->get_eta()) > 1.1)
	  { continue; }

	if(Verbosity() > 2)
	  {
	    std::cout << "Track with p= " << track->get_p() <<", eta / phi = "
		      << track->get_eta() << " / " << track->get_phi() 
		      << std::endl;
	  }

	_pflow_TRK_trk.push_back(track);
	_pflow_TRK_p.push_back(track->get_p());
	_pflow_TRK_eta.push_back(track->get_eta());
	_pflow_TRK_phi.push_back(track->get_phi());
	_pflow_TRK_match_EM.push_back( std::vector<int>() );
	_pflow_TRK_match_HAD.push_back( std::vector<int>() );
	_pflow_TRK_addtl_match_EM.push_back( std::vector< std::pair<int,float> >() );
	
	SvtxTrackState* cemcstate = track->get_state(cemcradius);
	SvtxTrackState* ohstate = track->get_state(ohcalradius);
	/// Get the track projections. If they failed for some reason, just use the track
	/// phi and eta values at the point of closest approach
	if(cemcstate)
	  {
	    _pflow_TRK_EMproj_phi.push_back(atan2(cemcstate->get_y(), cemcstate->get_x()));
	    _pflow_TRK_EMproj_eta.push_back(asinh(cemcstate->get_z()/sqrt(cemcstate->get_x()*cemcstate->get_x() + cemcstate->get_y()*cemcstate->get_y())));
	  } 
	  else { 
	  _pflow_TRK_EMproj_phi.push_back(track->get_phi());
	  _pflow_TRK_EMproj_eta.push_back(track->get_eta());
	  	  }
	if(ohstate)
	  {
	    _pflow_TRK_HADproj_phi.push_back(atan2(ohstate->get_py(), ohstate->get_px()));
	    _pflow_TRK_HADproj_eta.push_back(asinh(ohstate->get_pz()/ohstate->get_pt()));
	  }
	  else { 
	  _pflow_TRK_HADproj_phi.push_back(track->get_phi());
	  _pflow_TRK_HADproj_eta.push_back(track->get_eta());
	  }
      }

  } // 

  

  // read in EMCal topoClusters with E > 0.2 GeV
  {
    RawClusterContainer::ConstRange begin_end = clustersEM->getClusters();
    for ( RawClusterContainer::ConstIterator hiter = begin_end.first; hiter != begin_end.second; ++hiter)
      {
	float cluster_E = hiter->second->get_energy();
	if ( cluster_E < 0.2 ) continue;
	
	float cluster_phi = hiter->second->get_phi();
	/// default assume at vx_z = 0
	float cluster_theta = M_PI / 2.0 - atan2( hiter->second->get_z() , hiter->second->get_r() );
	float cluster_eta = -1 * log( tan( cluster_theta / 2.0 ) );
	
	if(vertex)
	  {
	    cluster_eta = RawClusterUtility::GetPseudorapidity(*(hiter->second),CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
	  }

	_pflow_EM_E.push_back( cluster_E );
	_pflow_EM_eta.push_back( cluster_eta );
	_pflow_EM_phi.push_back( cluster_phi );
	_pflow_EM_cluster.push_back(hiter->second);
	_pflow_EM_match_HAD.push_back( std::vector<int>() );
	_pflow_EM_match_TRK.push_back( std::vector<int>() );
	
	if ( Verbosity() > 5 && cluster_E > 0.2 ) 
	  std::cout << " EM topoCluster with E = " << cluster_E << ", eta / phi = " << cluster_eta << " / " << cluster_phi << " , nTow = " << hiter->second->getNTowers()  << std::endl;
	
	std::vector<float> this_cluster_tower_eta;
	std::vector<float> this_cluster_tower_phi;
	
	// read in towers
	RawCluster::TowerConstRange begin_end_towers = hiter->second->get_towers();
	for (RawCluster::TowerConstIterator iter = begin_end_towers.first; iter != begin_end_towers.second; ++iter) {
	  
	  if ( RawTowerDefs::decode_caloid( iter->first ) == RawTowerDefs::CalorimeterId::CEMC ) {
	    RawTower* tower = towersEM->getTower(iter->first);
	    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());
	    
	    this_cluster_tower_phi.push_back( tower_geom->get_phi() );
	    this_cluster_tower_eta.push_back( tower_geom->get_eta() );
	  }
	  else {
	    std::cout << "ParticleFlowReco::process_event : FATAL ERROR , EM topoClusters seem to contain HCal towers" << std::endl;
	    return Fun4AllReturnCodes::ABORTEVENT;
	  }	  
	} // close tower loop
	
	_pflow_EM_tower_eta.push_back( this_cluster_tower_eta );
	_pflow_EM_tower_phi.push_back( this_cluster_tower_phi );
	
      } // close cluster loop
    
  } // close 

  // read in HCal topoClusters with E > 0.2 GeV
  {
    RawClusterContainer::ConstRange begin_end = clustersHAD->getClusters();
    for ( RawClusterContainer::ConstIterator hiter = begin_end.first; hiter != begin_end.second; ++hiter)
      {
	float cluster_E = hiter->second->get_energy();
	if ( cluster_E < 0.2 ) continue;
	
	float cluster_phi = hiter->second->get_phi();
	// for now, assume event at vx_z = 0
	float cluster_theta = M_PI / 2.0 - atan2( hiter->second->get_z() , hiter->second->get_r() );
	float cluster_eta = -1 * log( tan( cluster_theta / 2.0 ) );
	if(vertex)
	  {
	    cluster_eta = RawClusterUtility::GetPseudorapidity(*(hiter->second),CLHEP::Hep3Vector(vertex->get_x(), vertex->get_y(), vertex->get_z()));
	  }

	_pflow_HAD_E.push_back( cluster_E );
	_pflow_HAD_eta.push_back( cluster_eta );
	_pflow_HAD_phi.push_back( cluster_phi );
	_pflow_HAD_cluster.push_back(hiter->second);

	_pflow_HAD_match_EM.push_back( std::vector<int>() );
	_pflow_HAD_match_TRK.push_back( std::vector<int>() );
	
	if ( Verbosity() > 5 && cluster_E > 0.2 ) 
	  std::cout << " HAD topoCluster with E = " << cluster_E << ", eta / phi = " << cluster_eta << " / " << cluster_phi << " , nTow = " << hiter->second->getNTowers()  << std::endl;
	
	std::vector<float> this_cluster_tower_eta;
	std::vector<float> this_cluster_tower_phi;
	
	// read in towers
	RawCluster::TowerConstRange begin_end_towers = hiter->second->get_towers();
	for (RawCluster::TowerConstIterator iter = begin_end_towers.first; iter != begin_end_towers.second; ++iter) {
	  
	  if ( RawTowerDefs::decode_caloid( iter->first ) == RawTowerDefs::CalorimeterId::HCALIN ) {

	    RawTower* tower = towersIH->getTower(iter->first);
	    RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
	    
	    this_cluster_tower_phi.push_back( tower_geom->get_phi() );
	    this_cluster_tower_eta.push_back( tower_geom->get_eta() );
	  }

	  else if ( RawTowerDefs::decode_caloid( iter->first ) == RawTowerDefs::CalorimeterId::HCALOUT ) {

	    RawTower* tower = towersOH->getTower(iter->first);
	    RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());
	    
	    this_cluster_tower_phi.push_back( tower_geom->get_phi() );
	    this_cluster_tower_eta.push_back( tower_geom->get_eta() );
	  } else {
	    std::cout << "ParticleFlowReco::process_event : FATAL ERROR , HCal topoClusters seem to contain EM towers" << std::endl;
	    return Fun4AllReturnCodes::ABORTEVENT;
	  }
	  
	  
	} // close tower loop
	
	_pflow_HAD_tower_eta.push_back( this_cluster_tower_eta );
	_pflow_HAD_tower_phi.push_back( this_cluster_tower_phi );
	
      } // close cluster loop
    
  } // close 

  // BEGIN LINKING STEP

  // Link TRK -> EM (best match, but keep reserve of others), and TRK -> HAD (best match)
  if ( Verbosity() > 2 ) 
    std::cout << "ParticleFlowReco::process_event : TRK -> EM and TRK -> HAD linking " << std::endl;

  for (unsigned int trk = 0; trk < _pflow_TRK_p.size() ; trk++ ) {
    
    if ( Verbosity() > 10 )
      std::cout << " TRK " << trk << " with p / eta / phi = " << _pflow_TRK_p[ trk ] << " / " << _pflow_TRK_eta[ trk ] << " / " << _pflow_TRK_phi[ trk ] << std::endl;

    // TRK -> EM link
    float min_em_dR = 0.2;
    int min_em_index = -1;
    
    for (unsigned int em = 0 ; em < _pflow_EM_E.size() ; em++) {

      float dR = calculate_dR( _pflow_TRK_EMproj_eta[ trk ] , _pflow_EM_eta[ em ] , _pflow_TRK_EMproj_phi[ trk ] , _pflow_EM_phi[ em ] );
    
      if ( dR > 0.2 ) continue;
   
      bool has_overlap = false;

      for (unsigned int tow = 0; tow < _pflow_EM_tower_eta.at( em ).size() ; tow++) {

	float tower_eta =  _pflow_EM_tower_eta.at( em ).at( tow );
	float tower_phi =  _pflow_EM_tower_phi.at( em ).at( tow );

	float deta = tower_eta - _pflow_TRK_EMproj_eta[ trk ];
	float dphi = tower_phi - _pflow_TRK_EMproj_phi[ trk ];
	if ( dphi > M_PI ) dphi -= 2 * M_PI;
	if ( dphi < -M_PI ) dphi += 2 * M_PI;

	if ( fabs( deta ) < 0.025 * 2.5 && fabs( dphi ) < 0.025 * 2.5 ) {
	  has_overlap = true;
	  break;
	}

      }

      if ( has_overlap ) {

	if ( Verbosity() > 5 ) 
	  std::cout << " -> possible match to EM " << em << " with dR = " << dR << std::endl;

	_pflow_TRK_addtl_match_EM.at( trk ).push_back( std::pair<int,float>( em, dR ) );

      } else {
	
	if ( Verbosity() > 5 ) 
	  std::cout << " -> no match to EM " << em << " (even though dR = " << dR << " )" << std::endl;
	
      }

    }

    // sort possible matches 

    std::sort( _pflow_TRK_addtl_match_EM.at( trk ).begin(), _pflow_TRK_addtl_match_EM.at( trk ).end(), sort_by_pair_second_lowest );
    if ( Verbosity() > 10 ) {
      for (unsigned int n = 0; n < _pflow_TRK_addtl_match_EM.at( trk ).size(); n++) {
	std::cout << " -> sorted list of matches, EM / dR = " <<  _pflow_TRK_addtl_match_EM.at( trk ).at( n ).first << " / " << _pflow_TRK_addtl_match_EM.at( trk ).at( n ).second << std::endl;
      }
    }
  
    if ( _pflow_TRK_addtl_match_EM.at( trk ).size() > 0 ) {
      min_em_index = _pflow_TRK_addtl_match_EM.at( trk ).at( 0 ).first;
      min_em_dR =  _pflow_TRK_addtl_match_EM.at( trk ).at( 0 ).second;
      // delete best matched element
      _pflow_TRK_addtl_match_EM.at( trk ).erase( _pflow_TRK_addtl_match_EM.at( trk ).begin() );
    }


    if ( min_em_index > -1 ) {
      _pflow_EM_match_TRK.at( min_em_index ).push_back( trk );
      _pflow_TRK_match_EM.at( trk ).push_back( min_em_index );

      if ( Verbosity() > 5 ) {
	std::cout << " -> matched EM " << min_em_index << " with pt / eta / phi = " << _pflow_EM_E.at( min_em_index ) << " / " << _pflow_EM_eta.at( min_em_index ) << " / " << _pflow_EM_phi.at( min_em_index ) << ", dR = " << min_em_dR;
	std::cout << " ( " << _pflow_TRK_addtl_match_EM.at( trk ).size() << " other possible matches ) " << std::endl;
      }
      
    } else {

      if ( Verbosity() > 5 ) 
	std::cout << " -> no EM match! ( best dR = " << min_em_dR << " ) " << std::endl;
    }
    
    // TRK -> HAD link
    float min_had_dR = 0.2;
    int min_had_index = -1;
    float max_had_pt = 0;

    // TODO: sequential linking should better happen here -- i.e. allow EM-matched HAD's into the possible pool
    for (unsigned int had = 0 ; had < _pflow_HAD_E.size() ; had++) {

      float dR = calculate_dR( _pflow_TRK_HADproj_eta[ trk ] , _pflow_HAD_eta[ had ] , _pflow_TRK_HADproj_phi[ trk ] , _pflow_HAD_phi[ had ] );
      
      if ( dR > 0.5 ) continue;
    
      bool has_overlap = false;

      for (unsigned int tow = 0; tow < _pflow_HAD_tower_eta.at( had ).size() ; tow++) {

	float tower_eta =  _pflow_HAD_tower_eta.at( had ).at( tow );
	float tower_phi =  _pflow_HAD_tower_phi.at( had ).at( tow );

	float deta = tower_eta - _pflow_TRK_HADproj_eta[ trk ];
	float dphi = tower_phi - _pflow_TRK_HADproj_phi[ trk ];
	if ( dphi > M_PI ) dphi -= 2 * M_PI;
	if ( dphi < -M_PI ) dphi += 2 * M_PI;

	if ( fabs( deta ) < 0.1 * 1.5 && fabs( dphi ) < 0.1 * 1.5 ) {
	  has_overlap = true;
	  break;
	}

      }

      if ( has_overlap ) {

	if ( Verbosity() > 5 ) 
	  std::cout << " -> possible match to HAD " << had << " with dR = " << dR << std::endl;

	if ( _pflow_HAD_E.at( had ) > max_had_pt ) {
	  max_had_pt = _pflow_HAD_E.at( had );
	  min_had_index = had;
	  min_had_dR = dR;
	}

      } else {
	
	if ( Verbosity() > 5 ) 
	  std::cout << " -> no match to HAD " << had << " (even though dR = " << dR << " )" << std::endl;
	
      }

    }

    if ( min_had_index > -1 ) {
      _pflow_HAD_match_TRK.at( min_had_index ).push_back( trk );
      _pflow_TRK_match_HAD.at( trk ).push_back( min_had_index );

      if ( Verbosity() > 5 ) 
	std::cout << " -> matched HAD " << min_had_index << " with pt / eta / phi = " << _pflow_HAD_E.at( min_had_index ) << " / " << _pflow_HAD_eta.at( min_had_index ) << " / " << _pflow_HAD_phi.at( min_had_index ) << ", dR = " << min_had_dR << std::endl;
      
    } else {
      if ( Verbosity() > 5 ) 
	std::cout << " -> no HAD match! ( best dR = " << min_had_dR << " ) " << std::endl;
    }


  }


  // EM->HAD linking 
  if ( Verbosity() > 2 ) 
    std::cout << "ParticleFlowReco::process_event : EM -> HAD linking " << std::endl;

  for (unsigned int em = 0; em < _pflow_EM_E.size() ; em++ ) {
    
    if ( Verbosity() > 10 )
      std::cout << " EM with E / eta / phi = " << _pflow_EM_E[ em ] << " / " << _pflow_EM_eta[ em ] << " / " << _pflow_EM_phi[ em ] << std::endl;
    
    // TRK -> HAD link
    float min_had_dR = 0.2;
    int min_had_index = -1;
    float max_had_pt = 0;
    
    for (unsigned int had = 0 ; had < _pflow_HAD_E.size() ; had++) {

      float dR = calculate_dR( _pflow_EM_eta[ em ] , _pflow_HAD_eta[ had ] , _pflow_EM_phi[ em ] , _pflow_HAD_phi[ had ] );
      if ( dR > 0.5 ) continue;

      bool has_overlap = false;

      for (unsigned int tow = 0; tow < _pflow_HAD_tower_eta.at( had ).size() ; tow++) {

	float tower_eta =  _pflow_HAD_tower_eta.at( had ).at( tow );
	float tower_phi =  _pflow_HAD_tower_phi.at( had ).at( tow );

	float deta = tower_eta - _pflow_EM_eta[ em ];
	float dphi = tower_phi - _pflow_EM_phi[ em ];
	if ( dphi > M_PI ) dphi -= 2 * M_PI;
	if ( dphi < -M_PI ) dphi += 2 * M_PI;

	if ( fabs( deta ) < 0.1 * 1.5 && fabs( dphi ) < 0.1 * 1.5 ) {
	  has_overlap = true;
	  break;
	}

      }

      if ( has_overlap ) {

	if ( Verbosity() > 5 ) 
	  std::cout << " -> possible match to HAD " << had << " with dR = " << dR << std::endl;

	if ( _pflow_HAD_E.at( had ) > max_had_pt ) {
	  max_had_pt = _pflow_HAD_E.at( had );
	  min_had_index = had;
	  min_had_dR = dR;
	}

      } else {
	
	if ( Verbosity() > 5 ) 
	  std::cout << " -> no match to HAD " << had << " (even though dR = " << dR << " )" << std::endl;
	
      }

    }

    if ( min_had_index > -1 ) {
      _pflow_HAD_match_EM.at( min_had_index ).push_back( em );
      _pflow_EM_match_HAD.at( em ).push_back( min_had_index );

      if ( Verbosity() > 5 ) 
	std::cout << " -> matched HAD with E / eta / phi = " << _pflow_HAD_E.at( min_had_index ) << " / " << _pflow_HAD_eta.at( min_had_index ) << " / " << _pflow_HAD_phi.at( min_had_index ) << ", dR = " << min_had_dR << std::endl;
      
    } else {
      if ( Verbosity() > 5 ) 
	std::cout << " -> no HAD match! ( best dR = " << min_had_dR << " ) " << std::endl;
    }

  }


  // SEQUENTIAL MATCHING: if TRK -> EM and EM -> HAD, ensure that TRK -> HAD 
  if ( Verbosity() > 2 ) 
    std::cout << "ParticleFlowReco::process_event : sequential TRK -> EM && EM -> HAD ==> TRK -> HAD matching " << std::endl;

  for (unsigned int trk = 0; trk < _pflow_TRK_p.size() ; trk++ ) {

    // go through all matched EMs
    for (unsigned int i = 0; i < _pflow_TRK_match_EM.at( trk ).size(); i++) {

      int em = _pflow_TRK_match_EM.at( trk ).at( i );

      // if this EM has a matched HAD...
      for (unsigned int j = 0; j < _pflow_EM_match_HAD.at( em ).size() ; j++) {

	int had = _pflow_EM_match_HAD.at( em ).at( j );

	// and the TRK is NOT matched to this HAD...
	bool is_trk_matched_to_HAD = false;
	for (unsigned int k = 0; k < _pflow_TRK_match_HAD.at( trk ).size(); k++) {
	  int existing_had = _pflow_TRK_match_HAD.at( trk ).at( k );
	  if ( had == existing_had ) is_trk_matched_to_HAD = true;
	}	
	
	// if this is the case, create TRK->HAD link
	if ( ! is_trk_matched_to_HAD ) {
	  _pflow_TRK_match_HAD.at( trk ).push_back( had );
	  _pflow_HAD_match_TRK.at( had ).push_back( trk );

	  if ( Verbosity() > 5 )  {
	    std::cout << " TRK " << trk << " with pt / eta / phi = " << _pflow_TRK_p.at( trk ) << " / " << _pflow_TRK_eta.at( trk ) << " / " << _pflow_TRK_phi.at( trk ) << std::endl;
	    std::cout << " -> sequential match to HAD " << had << " through EM " << j << std::endl;
	  }

	}

      } // close the HAD loop

    } // close the EM loop

  } // close the TRK loop

  // TRK->EM->HAD removal
  if ( Verbosity() > 2 ) 
    std::cout << "ParticleFlowReco::process_event : resolve TRK(s) + EM(s) -> HAD systems " << std::endl;
  
  for (unsigned int had = 0; had < _pflow_HAD_E.size() ; had++ ) {

    // only consider HAD with matched tracks ... others we will deal with later 
    if ( _pflow_HAD_match_TRK.at( had ).size() == 0 ) continue;

    if ( Verbosity() > 5 ) {
      std::cout << " HAD " << had << " with E / eta / phi = " << _pflow_HAD_E.at( had ) << " / " << _pflow_HAD_eta.at( had ) << " / " << _pflow_HAD_phi.at( had ) << std::endl;
    }
    
    // setup for Sum-pT^trk -> calo prediction
    float total_TRK_p = 0;
    float total_expected_E = 0;
    float total_expected_E_var = 0;

    // begin with this HAD calo energy
    float total_EMHAD_E = _pflow_HAD_E.at( had );

    std::vector<RawCluster*> matchedEClusters;

    // iterate over the EMs matched to this HAD 
    for (unsigned int j = 0; j < _pflow_HAD_match_EM.at( had ).size() ; j++ ) {

      int em = _pflow_HAD_match_EM.at( had ).at( j );

      // ensure there is at least one track matched to this EM
      if ( _pflow_EM_match_TRK.at( em ).size() == 0 ) continue;

      // add it to the total calo E 
      total_EMHAD_E += _pflow_EM_E.at( em );
      matchedEClusters.push_back(_pflow_EM_cluster.at(em));
      if ( Verbosity() > 5 ) {
	std::cout << " -> -> LINKED EM " << em << " with E / eta / phi = " <<  _pflow_EM_E.at( em ) << " / " <<  _pflow_EM_eta.at( em ) << " / " <<  _pflow_EM_phi.at( em ) << std::endl;
      }

    }
  
    // iterate over the TRKs matched to this HAD
    for (unsigned int j = 0 ; j <  _pflow_HAD_match_TRK.at( had ).size() ; j++) {

      int trk = _pflow_HAD_match_TRK.at( had ).at( j );

      if ( Verbosity() > 5 ) {
	std::cout << " -> -> LINKED TRK " << trk << " with p / eta / phi = " << _pflow_TRK_p.at( trk ) << " / " << _pflow_TRK_eta.at( trk ) << " / " << _pflow_TRK_phi.at( trk ) << std::endl;
      }
      
      total_TRK_p += _pflow_TRK_p.at( trk );

      std::pair<float, float> expected_signature = get_expected_signature( trk );

      float expected_E_mean = expected_signature.first;
      float expected_E_sigma = expected_signature.second;

      if ( Verbosity() > 5 ) {
	std::cout << " -> -> -> expected calo signature is " << expected_E_mean << " +/- " << expected_E_sigma << std::endl;
      }
      
      total_expected_E += expected_E_mean;
      total_expected_E_var += pow( expected_E_sigma , 2 );

      // add PFlow element for each track
      ParticleFlowElement *pflow = new ParticleFlowElementv1();
      
      // assume pion mass
      TLorentzVector tlv; 
      tlv.SetPtEtaPhiM( _pflow_TRK_p[ trk ] / cosh( _pflow_TRK_eta[ trk ] ) , _pflow_TRK_eta[ trk ] , _pflow_TRK_phi[ trk ] , 0.135 ); 

      pflow->set_px( tlv.Px() );
      pflow->set_py( tlv.Py() );
      pflow->set_pz( tlv.Pz() );
      pflow->set_e( tlv.E() );
      pflow->set_track(_pflow_TRK_trk[ trk ]);
      pflow->set_eclusters(matchedEClusters);
      pflow->set_hcluster(_pflow_HAD_cluster.at(had));
      pflow->set_id( global_pflow_index );
      pflow->set_type( ParticleFlowElement::PFLOWTYPE::MATCHED_CHARGED_HADRON );
    
      pflowContainer->AddParticleFlowElement( global_pflow_index, pflow );
      global_pflow_index++;

    }
    // Track + E+HCal PF elements are created
    
    // process compatibility of fit 
    float total_expected_E_err = sqrt( total_expected_E_var );

    if ( Verbosity() > 5 ) {
      std::cout << " -> Total track Sum p = " << total_TRK_p << " , expected calo Sum E = " << total_expected_E << " +/- " << total_expected_E_err << " , observed EM+HAD Sum E = " << total_EMHAD_E << std::endl;
    }

    // if Sum pT > calo, add in additional possible matched EMs associated with tracks until that is no longer the case

    if ( total_expected_E > total_EMHAD_E ) {
      
      if ( Verbosity() > 5 ) 
	std::cout << " -> Expected E > Observed E, looking for additional potential TRK->EM matches" << std::endl;
      
      std::map<int, float> additional_EMs;
      
      for (unsigned int j = 0 ; j <  _pflow_HAD_match_TRK.at( had ).size() ; j++) {
	
	int trk = _pflow_HAD_match_TRK.at( had ).at( j );

	int addtl_matches = _pflow_TRK_addtl_match_EM.at( trk ).size();

	if ( Verbosity() > 10 )
	  std::cout << " -> -> TRK " << trk << " has " << addtl_matches << " additional matches! " << std::endl;
	
	for (unsigned int n  = 0 ; n < _pflow_TRK_addtl_match_EM.at( trk ).size() ; n++ ) {
	  if ( Verbosity() > 10 )
	    std::cout << " -> -> -> additional match to EM = " << _pflow_TRK_addtl_match_EM.at( trk ).at( n ).first << " with dR = " <<  _pflow_TRK_addtl_match_EM.at( trk ).at( n ).second << std::endl;
	  
	  float existing_dR = 0.21;
	  int counts = additional_EMs.count(  _pflow_TRK_addtl_match_EM.at( trk ).at( n ).first );
	  if ( counts > 0 ) {
	    existing_dR = additional_EMs[ _pflow_TRK_addtl_match_EM.at( trk ).at( n ).first ];
	  }
	  if ( _pflow_TRK_addtl_match_EM.at( trk ).at( n ).second < existing_dR )
	    additional_EMs[ _pflow_TRK_addtl_match_EM.at( trk ).at( n ).first ] = _pflow_TRK_addtl_match_EM.at( trk ).at( n ).second;
	}
	
      }

      // map now assured to have only minimal dR values for each possible additional EM
      // translate the map to a vector of pairs, then sort by smallest dR

      std::vector< std::pair<int,float> > additional_EMs_vec;

      for (auto& x : additional_EMs) {
	additional_EMs_vec.push_back( std::pair<int,float>( x.first , x.second ) );
      }

      std::sort( additional_EMs_vec.begin(), additional_EMs_vec.end(), sort_by_pair_second_lowest );

      if ( Verbosity() > 5 )      
	std::cout << " -> Sorting the set of potential additional EMs " << std::endl;

      // now add in additional EMs until there are none left or it is no longer the case that Sum pT > calo

      int n_EM_added = 0;
      while ( additional_EMs_vec.size() != 0 && total_expected_E > total_EMHAD_E ) {

	int new_EM = additional_EMs_vec.at( 0 ).first;

	if ( Verbosity() > 5 )      
	  std::cout << " -> adding EM " << new_EM << " ( dR = " << additional_EMs_vec.at( 0 ).second << " to the system (should not see it as orphan below)" << std::endl;

	// for now, just make the first HAD-linked track point to this new EM, and vice versa
	_pflow_EM_match_TRK.at( new_EM ).push_back( _pflow_HAD_match_TRK.at( had ).at( 0 ) );
	_pflow_TRK_match_EM.at( _pflow_HAD_match_TRK.at( had ).at( 0 ) ).push_back( new_EM );
	
	// add to expected calo
	total_EMHAD_E += _pflow_EM_E.at( new_EM );

	// erase lowest-dR EM
	additional_EMs_vec.erase( additional_EMs_vec.begin() );

	n_EM_added++;
      }
    
      if ( Verbosity() > 5) {
	if ( n_EM_added > 0 ) {
	  std::cout << "After adding N = " << n_EM_added << " any additional EMs : " << std::endl;
	  std::cout << "-> Total track Sum p = " << total_TRK_p << " , expected calo Sum E = " << total_expected_E << " +/- " << total_expected_E_err << " , observed EM+HAD Sum E = " << total_EMHAD_E << std::endl;
	}
	else { 
	  std::cout << "No additional EMs found, continuing hypothesis check" << std::endl;
	}
      }
    }
    
  
  
    if ( total_expected_E + _energy_match_Nsigma * total_expected_E_err > total_EMHAD_E ) {
      
      if ( Verbosity() > 5 ) {
	std::cout << " -> -> calo compatible within Nsigma = " << _energy_match_Nsigma << " , remove and keep tracks " << std::endl;
      }

      // PFlow elements already created from tracks above, no more needs to be done
      
    } else {

      float residual_energy = total_EMHAD_E - total_expected_E;

      if ( Verbosity() > 5 ) {
	std::cout << " -> -> calo not compatible, create leftover cluster with " << residual_energy << std::endl;
      }

      // create additional PFlow element (tracks already created above)
      ParticleFlowElement *pflow = new ParticleFlowElementv1();
      
      // assume no mass, but could update to use K0L mass(?)
      TLorentzVector tlv; tlv.SetPtEtaPhiM( residual_energy / cosh( _pflow_HAD_eta[ had ] ) , _pflow_HAD_eta[ had ] , _pflow_HAD_phi[ had ] , 0 ); 

      pflow->set_px( tlv.Px() );
      pflow->set_py( tlv.Py() );
      pflow->set_pz( tlv.Pz() );
      pflow->set_e( tlv.E() );
      pflow->set_track(nullptr);
      pflow->set_eclusters(matchedEClusters);
      pflow->set_hcluster(_pflow_HAD_cluster.at(had));
      pflow->set_id( global_pflow_index );
      pflow->set_type( ParticleFlowElement::PFLOWTYPE::LEFTOVER_EM_PARTICLE );

      pflowContainer->AddParticleFlowElement( global_pflow_index, pflow );
      global_pflow_index++;

    }

  } // close HAD loop 

  // TRK->EM removal

  if ( Verbosity() > 2 ) 
    std::cout << "ParticleFlowReco::process_event : resolve TRK(s) -> EM(s) ( + no HAD) systems " << std::endl;
  
  for (unsigned int em = 0; em < _pflow_EM_E.size() ; em++ ) {

    // only consider EM with matched tracks, but no matched HADs
    if ( _pflow_EM_match_HAD.at( em ).size() != 0 ) continue;
    if ( _pflow_EM_match_TRK.at( em ).size() == 0 ) continue;

    if ( Verbosity() > 5 ) {
      std::cout << " EM " << em << " with E / eta / phi = " << _pflow_EM_E.at( em ) << " / " << _pflow_EM_eta.at( em ) << " / " << _pflow_EM_phi.at( em ) << std::endl;
    }
    
    // setup for Sum-pT^trk -> calo prediction
    float total_TRK_p = 0;
    float total_expected_E = 0;
    float total_expected_E_var = 0;

    // begin with this EM calo energy
    float total_EM_E = _pflow_EM_E.at( em );
  
    // iterate over the TRKs matched to this EM
    for (unsigned int j = 0 ; j <  _pflow_EM_match_TRK.at( em ).size() ; j++) {

      int trk = _pflow_EM_match_TRK.at( em ).at( j );

      if ( Verbosity() > 5 ) {
	std::cout << " -> -> LINKED TRK with p / eta / phi = " << _pflow_TRK_p.at( trk ) << " / " << _pflow_TRK_eta.at( trk ) << " / " << _pflow_TRK_phi.at( trk ) << std::endl;
      }
      
      total_TRK_p += _pflow_TRK_p.at( trk );

      std::pair<float, float> expected_signature = get_expected_signature( trk );

      float expected_E_mean = expected_signature.first;
      float expected_E_sigma = expected_signature.second;

      if ( Verbosity() > 5 ) {
	std::cout << " -> -> -> expected calo signature is " << expected_E_mean << " +/- " << expected_E_sigma << std::endl;
      }
      
      total_expected_E += expected_E_mean;
      total_expected_E_var += pow( expected_E_sigma , 2 );

      // add PFlow element for each track
      ParticleFlowElement *pflow = new ParticleFlowElementv1();
      
      // assume pion mass
      TLorentzVector tlv; tlv.SetPtEtaPhiM( _pflow_TRK_p[ trk ] / cosh( _pflow_TRK_eta[ trk ] ) , _pflow_TRK_eta[ trk ] , _pflow_TRK_phi[ trk ] , 0.135 ); 

      std::vector<RawCluster*> eclus;
      eclus.push_back(_pflow_EM_cluster.at(em));

      pflow->set_px( tlv.Px() );
      pflow->set_py( tlv.Py() );
      pflow->set_pz( tlv.Pz() );
      pflow->set_e( tlv.E() );
      pflow->set_track(_pflow_TRK_trk.at(trk));
      pflow->set_eclusters(eclus);
      pflow->set_hcluster(nullptr);
      pflow->set_id( global_pflow_index );
      pflow->set_type( ParticleFlowElement::PFLOWTYPE::MATCHED_CHARGED_HADRON );

      pflowContainer->AddParticleFlowElement( global_pflow_index, pflow );
      global_pflow_index++;

    }

    // process compatibility of fit 
    float total_expected_E_err = sqrt( total_expected_E_var );

    if ( Verbosity() > 5 ) {
      std::cout << " -> Total track Sum p = " << total_TRK_p << " , expected calo Sum E = " << total_expected_E << " +/- " << total_expected_E_err << " , observed EM Sum E = " << total_EM_E << std::endl;
    }
      
    if ( total_expected_E + _energy_match_Nsigma * total_expected_E_err > total_EM_E ) {
      
      if ( Verbosity() > 5 ) {
	std::cout << " -> -> calo compatible within Nsigma = " << _energy_match_Nsigma << "  , remove and keep tracks " << std::endl;
      }

      // PFlow elements already created from tracks above, no more needs to be done
      
    } else {

      float residual_energy = total_EM_E - total_expected_E;

      if ( Verbosity() > 5 ) {
	std::cout << " -> -> calo not compatible, create leftover cluster with " << residual_energy << std::endl;
      }

      // create additional PFlow element (tracks already created above)
      ParticleFlowElement *pflow = new ParticleFlowElementv1();
      
      // assume no mass, but could update to use K0L mass(?)
      TLorentzVector tlv; tlv.SetPtEtaPhiM( residual_energy / cosh( _pflow_EM_eta[ em ] ) , _pflow_EM_eta[ em ] , _pflow_EM_phi[ em ] , 0 ); 

      std::vector<RawCluster*> eclus;
      eclus.push_back(_pflow_EM_cluster.at(em));

      pflow->set_px( tlv.Px() );
      pflow->set_py( tlv.Py() );
      pflow->set_pz( tlv.Pz() );
      pflow->set_e( tlv.E() );
      pflow->set_eclusters(eclus);
      pflow->set_hcluster(nullptr);
      pflow->set_track(nullptr);
      pflow->set_id( global_pflow_index );
      pflow->set_type( ParticleFlowElement::PFLOWTYPE::LEFTOVER_EM_PARTICLE );

      pflowContainer->AddParticleFlowElement( global_pflow_index, pflow );
      global_pflow_index++;

    }

  } // close EM loop 

  // now remove unmatched elements

  if ( Verbosity() > 2 ) 
    std::cout << "ParticleFlowReco::process_event : remove TRK-unlinked EMs and HADs " << std::endl;
  
  for (unsigned int em = 0; em < _pflow_EM_E.size() ; em++ ) {

    // only consider EMs withOUT matched tracks ... we have dealt with the matched cases above
    if ( _pflow_EM_match_TRK.at( em ).size() != 0 ) continue;

    if ( Verbosity() > 5 ) {
      std::cout << " unmatched EM " << em << " with E / eta / phi = " << _pflow_EM_E.at( em ) << " / " << _pflow_EM_eta.at( em ) << " / " << _pflow_EM_phi.at( em ) << std::endl;
    }
    
    // add PFlow element for this EM 
    ParticleFlowElement *pflow = new ParticleFlowElementv1();
    
    // assume massless, could be updated to use K0L
    TLorentzVector tlv; tlv.SetPtEtaPhiM( _pflow_EM_E[ em ] / cosh( _pflow_EM_eta[ em ] ) , _pflow_EM_eta[ em ] , _pflow_EM_phi[ em ] , 0 ); 
    
    std::vector<RawCluster*> eclus;
    eclus.push_back(_pflow_EM_cluster.at(em));

    pflow->set_px( tlv.Px() );
    pflow->set_py( tlv.Py() );
    pflow->set_pz( tlv.Pz() );
    pflow->set_e( tlv.E() );
    pflow->set_eclusters(eclus);
    pflow->set_hcluster(nullptr);
    pflow->set_track(nullptr);
    pflow->set_id( global_pflow_index );
    pflow->set_type( ParticleFlowElement::PFLOWTYPE::UNMATCHED_EM_PARTICLE );

    pflowContainer->AddParticleFlowElement( global_pflow_index, pflow );
    global_pflow_index++;
    
    
  } // close EM loop 

  for (unsigned int had = 0; had < _pflow_HAD_E.size() ; had++ ) {

    // only consider HADs withOUT matched tracks ... we have dealt with the matched cases above
    if ( _pflow_HAD_match_TRK.at( had ).size() != 0 ) continue;

    if ( Verbosity() > 5 ) {
      std::cout << " unmatched HAD " << had << " with E / eta / phi = " << _pflow_HAD_E.at( had ) << " / " << _pflow_HAD_eta.at( had ) << " / " << _pflow_HAD_phi.at( had ) << std::endl;
    }
    
    // add PFlow element for this HAD 
    ParticleFlowElement *pflow = new ParticleFlowElementv1();
    
    // assume massless, could be updated to use K0L
    TLorentzVector tlv; tlv.SetPtEtaPhiM( _pflow_HAD_E[ had ] / cosh( _pflow_HAD_eta[ had ] ) , _pflow_HAD_eta[ had ] , _pflow_HAD_phi[ had ] , 0 ); 
    
    pflow->set_px( tlv.Px() );
    pflow->set_py( tlv.Py() );
    pflow->set_pz( tlv.Pz() );
    pflow->set_e( tlv.E() );
    pflow->set_track(nullptr);
    pflow->set_eclusters(std::vector<RawCluster*>());
    pflow->set_hcluster(_pflow_HAD_cluster.at(had));
    pflow->set_id( global_pflow_index );
    pflow->set_type( ParticleFlowElement::PFLOWTYPE::UNMATCHED_NEUTRAL_HADRON );

    pflowContainer->AddParticleFlowElement( global_pflow_index, pflow );
    global_pflow_index++;
    
    
  } // close HAD loop 

  for (unsigned int trk = 0; trk < _pflow_TRK_p.size() ; trk++ ) {

    // only consider TRKs withOUT matched EM or HAD 
    if ( _pflow_TRK_match_EM.at( trk ).size() != 0 || _pflow_TRK_match_HAD.at( trk ).size() != 0 ) continue;

    if ( Verbosity() > 5 ) {
      std::cout << " unmatched TRK " << trk << " with p / eta / phi = " << _pflow_TRK_p.at( trk ) << " / " << _pflow_TRK_eta.at( trk ) << " / " << _pflow_TRK_phi.at( trk ) << std::endl;
    }
    
    // add PFlow element for this TRK 
    ParticleFlowElement *pflow = new ParticleFlowElementv1();
    
    // assume massless, could be updated to use K0L
    TLorentzVector tlv; tlv.SetPtEtaPhiM( _pflow_TRK_p[ trk ] / cosh( _pflow_TRK_eta[ trk ] ) , _pflow_TRK_eta[ trk ] , _pflow_TRK_phi[ trk ] , 0.135 ); 
    
    pflow->set_px( tlv.Px() );
    pflow->set_py( tlv.Py() );
    pflow->set_pz( tlv.Pz() );
    pflow->set_e( tlv.E() );
    pflow->set_track(_pflow_TRK_trk.at(trk));
    pflow->set_eclusters(std::vector<RawCluster*>());
    pflow->set_hcluster(nullptr);
    pflow->set_id( global_pflow_index );
    pflow->set_type( ParticleFlowElement::PFLOWTYPE::UNMATCHED_CHARGED_HADRON );

    pflowContainer->AddParticleFlowElement( global_pflow_index, pflow );
    global_pflow_index++;
        
  } // close TRK loop 

  // DEBUG: print out all PFLow elements
  if ( Verbosity() > 5 ) {
    std::cout << "ParticleFlowReco::process_event: summary of PFlow objects " << std::endl;
    
    ParticleFlowElementContainer::ConstRange begin_end = pflowContainer->getParticleFlowElements();
    for ( ParticleFlowElementContainer::ConstIterator hiter = begin_end.first; hiter != begin_end.second; ++hiter) {
      hiter->second->identify();
    }

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int ParticleFlowReco::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // store the PFlow elements under a sub-node directory
  PHCompositeNode *pflowNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PARTICLEFLOW"));
  if (!pflowNode)
    {
      pflowNode = new PHCompositeNode("PARTICLEFLOW");
      dstNode->addNode( pflowNode );
    }

  // create the ParticleFlowElementContainer node...
  ParticleFlowElementContainer *pflowElementContainer = findNode::getClass<ParticleFlowElementContainer>(topNode, "ParticleFlowElements");
  if (!pflowElementContainer)
    {
      pflowElementContainer = new ParticleFlowElementContainer();
      PHIODataNode<PHObject> *pflowElementNode = new PHIODataNode<PHObject>(pflowElementContainer, "ParticleFlowElements", "PHObject");
      pflowNode->addNode( pflowElementNode );
    }
  else
    {
      std::cout << PHWHERE << "::ERROR - ParticleFlowElements node alerady exists, but should not" << std::endl;
      exit(-1);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ParticleFlowReco::ResetEvent(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 0)
  {
  std::cout << "ParticleFlowReco::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ParticleFlowReco::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
  {
  std::cout << "ParticleFlowReco::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ParticleFlowReco::End(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 0)
  {
  std::cout << "ParticleFlowReco::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ParticleFlowReco::Reset(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 0)
  {
 std::cout << "ParticleFlowReco::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void ParticleFlowReco::Print(const std::string &what) const
{
  std::cout << "ParticleFlowReco::Print(const std::string &what) const Printing info for " << what << std::endl;
}
