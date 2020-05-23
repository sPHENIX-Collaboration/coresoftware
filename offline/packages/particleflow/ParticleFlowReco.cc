#include "ParticleFlowReco.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include "TLorentzVector.h"
#include <iostream>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <phool/getClass.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include "ParticleFlowElementContainer.h"
#include "ParticleFlowElementv1.h"

float ParticleFlowReco::calculate_dR( float eta1, float eta2, float phi1, float phi2 ) {

  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  while ( dphi > 3.14159 ) dphi -= 2 * 3.14159;
  while ( dphi < -3.14159 ) dphi += 2 * 3.14159;
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
 SubsysReco(name)
{
  std::cout << "ParticleFlowReco::ParticleFlowReco(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
ParticleFlowReco::~ParticleFlowReco()
{
  std::cout << "ParticleFlowReco::~ParticleFlowReco() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int ParticleFlowReco::Init(PHCompositeNode *topNode)
{
  std::cout << "ParticleFlowReco::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ParticleFlowReco::InitRun(PHCompositeNode *topNode)
{
  std::cout << "ParticleFlowReco::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;

  return CreateNode(topNode);

}

//____________________________________________________________________________..
int ParticleFlowReco::process_event(PHCompositeNode *topNode)
{

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

  _pflow_EM_E.clear();
  _pflow_EM_eta.clear();
  _pflow_EM_phi.clear();
  _pflow_EM_tower_eta.clear();
  _pflow_EM_tower_phi.clear();
  _pflow_EM_match_HAD.clear();
  _pflow_EM_match_TRK.clear();

  _pflow_HAD_E.clear();
  _pflow_HAD_eta.clear();
  _pflow_HAD_phi.clear();
  _pflow_HAD_tower_eta.clear();
  _pflow_HAD_tower_phi.clear();
  _pflow_HAD_match_EM.clear();
  _pflow_HAD_match_TRK.clear();


  if ( Verbosity() > 2 ) 
    std::cout << "ParticleFlowReco::process_event : initial population of TRK, EM, HAD objects " << std::endl;

  // read in tracks with > 0.5 GeV
  // currently just do this from truth
  {
    
    PHG4TruthInfoContainer* truthinfo = findNode::getClass <PHG4TruthInfoContainer> (topNode, "G4TruthInfo");
    PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange ();
    
    for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter) {
      PHG4Particle* g4particle = iter->second;
      
      TLorentzVector t;
      t.SetPxPyPzE( g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e() );
      
      float truth_pt = t.Pt();
      float truth_p = t.P();
      if (truth_pt < 0.5)
	continue; // only keep pt > 0.5 GeV
      
      float truth_eta = t.Eta();
      if (fabs (truth_eta) > 1.1)
	continue; // only keep |eta| < 1.1
      
      float truth_phi = t.Phi();
      
      int truth_pid = abs( g4particle->get_pid() ); // particle species

      // only keep charged truth particles
      if ( truth_pid != 211 && truth_pid != 321 && truth_pid != 2212 ) continue;

      _pflow_TRK_p.push_back( truth_p );
      _pflow_TRK_eta.push_back( truth_eta );
      _pflow_TRK_phi.push_back( truth_phi );

      _pflow_TRK_match_EM.push_back( std::vector<int>() );
      _pflow_TRK_match_HAD.push_back( std::vector<int>() );

      if ( Verbosity() > 5 && truth_pt > 0.5 ) 
	std::cout << " TRK with p / pT = " << truth_p << " / " << truth_pt  << " , eta / phi = " << truth_eta << " / " << truth_phi << std::endl;
      
    } // close truth paticle loop

  } // 

  

  // read in EMCal topoClusters with E > 0.2 GeV
  {
    RawClusterContainer::ConstRange begin_end = clustersEM->getClusters();
    for ( RawClusterContainer::ConstIterator hiter = begin_end.first; hiter != begin_end.second; ++hiter)
      {
	float cluster_E = hiter->second->get_energy();
	if ( cluster_E < 0.2 ) continue;
	
	float cluster_phi = hiter->second->get_phi();
	// for now, assume event at vx_z = 0
	float cluster_theta = 3.14159 / 2.0 - atan2( hiter->second->get_z() , hiter->second->get_r() );
	float cluster_eta = -1 * log( tan( cluster_theta / 2.0 ) );
	
	_pflow_EM_E.push_back( cluster_E );
	_pflow_EM_eta.push_back( cluster_eta );
	_pflow_EM_phi.push_back( cluster_phi );
	
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
	float cluster_theta = 3.14159 / 2.0 - atan2( hiter->second->get_z() , hiter->second->get_r() );
	float cluster_eta = -1 * log( tan( cluster_theta / 2.0 ) );
	
	_pflow_HAD_E.push_back( cluster_E );
	_pflow_HAD_eta.push_back( cluster_eta );
	_pflow_HAD_phi.push_back( cluster_phi );
	
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
      std::cout << " TRK with p / eta / phi = " << _pflow_TRK_p[ trk ] << " / " << _pflow_TRK_eta[ trk ] << " / " << _pflow_TRK_phi[ trk ] << std::endl;

    // TRK -> EM link
    float min_em_dR = 0.2;
    int min_em_index = -1;
    float max_em_pt = 0;
    
    for (unsigned int em = 0 ; em < _pflow_EM_E.size() ; em++) {

      float dR = calculate_dR( _pflow_TRK_eta[ trk ] , _pflow_EM_eta[ em ] , _pflow_TRK_phi[ trk ] , _pflow_EM_phi[ em ] );
      if ( dR > 0.2 ) continue;

      bool has_overlap = false;

      for (unsigned int tow = 0; tow < _pflow_EM_tower_eta.at( em ).size() ; tow++) {

	float tower_eta =  _pflow_EM_tower_eta.at( em ).at( tow );
	float tower_phi =  _pflow_EM_tower_phi.at( em ).at( tow );

	float deta = tower_eta - _pflow_TRK_eta[ trk ];
	float dphi = tower_phi - _pflow_TRK_phi[ trk ];
	if ( dphi > 3.14159 ) dphi -= 2 * 3.14159;
	if ( dphi < -3.14159 ) dphi += 2 * 3.14159;

	if ( fabs( deta ) < 0.025 * 2.5 && fabs( dphi ) < 0.025 * 2.5 ) {
	  has_overlap = true;
	  break;
	}

      }

      if ( has_overlap ) {

	if ( Verbosity() > 5 ) 
	  std::cout << " -> possible match to EM " << em << " with dR = " << dR << std::endl;

	if ( _pflow_EM_E.at( em ) > max_em_pt ) {
	  max_em_pt = _pflow_EM_E.at( em );
	  min_em_index = em;
	  min_em_dR = dR;
	}

      } else {
	
	if ( Verbosity() > 5 ) 
	  std::cout << " -> no match to EM " << em << " (even though dR = " << dR << " )" << std::endl;
	
      }

    }

    if ( min_em_index > -1 ) {
      _pflow_EM_match_TRK.at( min_em_index ).push_back( trk );
      _pflow_TRK_match_EM.at( trk ).push_back( min_em_index );

      if ( Verbosity() > 5 ) 
	std::cout << " -> matched EM with pt / eta / phi = " << _pflow_EM_E.at( min_em_index ) << " / " << _pflow_EM_eta.at( min_em_index ) << " / " << _pflow_EM_phi.at( min_em_index ) << ", dR = " << min_em_dR << std::endl;
      
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

      float dR = calculate_dR( _pflow_TRK_eta[ trk ] , _pflow_HAD_eta[ had ] , _pflow_TRK_phi[ trk ] , _pflow_HAD_phi[ had ] );
      if ( dR > 0.5 ) continue;

      bool has_overlap = false;

      for (unsigned int tow = 0; tow < _pflow_HAD_tower_eta.at( had ).size() ; tow++) {

	float tower_eta =  _pflow_HAD_tower_eta.at( had ).at( tow );
	float tower_phi =  _pflow_HAD_tower_phi.at( had ).at( tow );

	float deta = tower_eta - _pflow_TRK_eta[ trk ];
	float dphi = tower_phi - _pflow_TRK_phi[ trk ];
	if ( dphi > 3.14159 ) dphi -= 2 * 3.14159;
	if ( dphi < -3.14159 ) dphi += 2 * 3.14159;

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
	std::cout << " -> matched HAD with pt / eta / phi = " << _pflow_HAD_E.at( min_had_index ) << " / " << _pflow_HAD_eta.at( min_had_index ) << " / " << _pflow_HAD_phi.at( min_had_index ) << ", dR = " << min_had_dR << std::endl;
      
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
	if ( dphi > 3.14159 ) dphi -= 2 * 3.14159;
	if ( dphi < -3.14159 ) dphi += 2 * 3.14159;

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

	  std::cout << " TRK " << trk << " with pt / eta / phi = " << _pflow_TRK_p.at( trk ) << " / " << _pflow_TRK_eta.at( trk ) << " / " << _pflow_TRK_phi.at( trk ) << std::endl;
	  std::cout << " -> sequential match to HAD " << had << " through EM " << j << std::endl;

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

    // iterate over the EMs matched to this HAD 
    for (unsigned int j = 0; j < _pflow_HAD_match_EM.at( had ).size() ; j++ ) {

      int em = _pflow_HAD_match_EM.at( had ).at( j );

      // ensure there is at least one track matched to this EM
      if ( _pflow_EM_match_TRK.at( em ).size() == 0 ) continue;

      // add it to the total calo E 
      total_EMHAD_E += _pflow_EM_E.at( em );
      
      if ( Verbosity() > 5 ) {
	std::cout << " -> -> LINKED EM " << em << " with E / eta / phi = " <<  _pflow_EM_E.at( em ) << " / " <<  _pflow_EM_eta.at( em ) << " / " <<  _pflow_EM_phi.at( em ) << std::endl;
      }

    }
  
    // iterate over the TRKs matched to this HAD
    for (unsigned int j = 0 ; j <  _pflow_HAD_match_TRK.at( had ).size() ; j++) {

      int trk = _pflow_HAD_match_TRK.at( had ).at( j );

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

      pflow->set_px( tlv.Px() );
      pflow->set_py( tlv.Py() );
      pflow->set_pz( tlv.Pz() );
      pflow->set_e( tlv.E() );
      pflow->set_id( global_pflow_index );

      pflowContainer->AddParticleFlowElement( global_pflow_index, pflow );
      global_pflow_index++;

    }

    // process compatibility of fit 
    float total_expected_E_err = sqrt( total_expected_E_var );

    if ( Verbosity() > 5 ) {
      std::cout << " -> Total track Sum p = " << total_TRK_p << " , expected calo Sum E = " << total_expected_E << " +/- " << total_expected_E_err << " , observed EM+HAD Sum E = " << total_EMHAD_E << std::endl;
    }
      
    if ( total_expected_E + 1.5 * total_expected_E_err > total_EMHAD_E ) {
      
      if ( Verbosity() > 5 ) {
	std::cout << " -> -> calo compatible within 1.5 sigma, remove and keep tracks " << std::endl;
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
      pflow->set_id( global_pflow_index );

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

      pflow->set_px( tlv.Px() );
      pflow->set_py( tlv.Py() );
      pflow->set_pz( tlv.Pz() );
      pflow->set_e( tlv.E() );
      pflow->set_id( global_pflow_index );

      pflowContainer->AddParticleFlowElement( global_pflow_index, pflow );
      global_pflow_index++;

    }

    // process compatibility of fit 
    float total_expected_E_err = sqrt( total_expected_E_var );

    if ( Verbosity() > 5 ) {
      std::cout << " -> Total track Sum p = " << total_TRK_p << " , expected calo Sum E = " << total_expected_E << " +/- " << total_expected_E_err << " , observed EM Sum E = " << total_EM_E << std::endl;
    }
      
    if ( total_expected_E + 1.5 * total_expected_E_err > total_EM_E ) {
      
      if ( Verbosity() > 5 ) {
	std::cout << " -> -> calo compatible within 1.5 sigma, remove and keep tracks " << std::endl;
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

      pflow->set_px( tlv.Px() );
      pflow->set_py( tlv.Py() );
      pflow->set_pz( tlv.Pz() );
      pflow->set_e( tlv.E() );
      pflow->set_id( global_pflow_index );

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
    
    pflow->set_px( tlv.Px() );
    pflow->set_py( tlv.Py() );
    pflow->set_pz( tlv.Pz() );
    pflow->set_e( tlv.E() );
    pflow->set_id( global_pflow_index );
    
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
    pflow->set_id( global_pflow_index );

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
    pflow->set_id( global_pflow_index );

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



  std::cout << "ParticleFlowReco::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
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
int ParticleFlowReco::ResetEvent(PHCompositeNode *topNode)
{
  std::cout << "ParticleFlowReco::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ParticleFlowReco::EndRun(const int runnumber)
{
  std::cout << "ParticleFlowReco::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ParticleFlowReco::End(PHCompositeNode *topNode)
{
  std::cout << "ParticleFlowReco::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ParticleFlowReco::Reset(PHCompositeNode *topNode)
{
 std::cout << "ParticleFlowReco::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void ParticleFlowReco::Print(const std::string &what) const
{
  std::cout << "ParticleFlowReco::Print(const std::string &what) const Printing info for " << what << std::endl;
}
