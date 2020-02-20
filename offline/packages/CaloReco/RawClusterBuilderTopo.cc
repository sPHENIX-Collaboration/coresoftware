#include "RawClusterBuilderTopo.h"

#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawClusterDefs.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/phool.h>

#include <cassert>
#include <cmath>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include <algorithm>

bool sort_by_pair_second( const std::pair<int,float> &a,  const std::pair<int,float> &b) 
{ 
  return (a.second > b.second); 
} 

int RawClusterBuilderTopo::RawClusterBuilderTopo_constants_EMCal_eta_start_given_IHCal[24] = 
  {2,  6,  10, 14, 18, 22, 26, 30, 33, 37, 41, 44, 
   48, 52, 55, 59, 63, 66, 70, 74, 78, 82, 86, 90 };

int RawClusterBuilderTopo::RawClusterBuilderTopo_constants_EMCal_eta_end_given_IHCal[24] = 
  {5,  9,  13, 17, 21, 25, 29, 32, 36, 40, 43, 47,
   51, 54, 58, 62, 65, 69, 73, 77, 81, 85, 89, 93 };

int RawClusterBuilderTopo::RawClusterBuilderTopo_constants_IHCal_eta_given_EMCal[96] = {
  -1, -1,  0,  0,  0,  0,  1,  1,  1,  1,  2,  2, 
   2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  5,  5, 
   5,  5,  6,  6,  6,  6,  7,  7,  7,  8,  8,  8, 
   8,  9,  9,  9,  9, 10, 10, 10, 11, 11, 11, 11, 
  12, 12, 12, 12, 13, 13, 13, 14, 14, 14, 14, 15, 
  15, 15, 15, 16, 16, 16, 17, 17, 17, 17, 18, 18, 
  18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 
  21, 21, 22, 22, 22, 22, 23, 23, 23, 23, -1, -1 };

std::vector<int> RawClusterBuilderTopo::get_adjacent_towers_by_ID( int ID ) {
  
  int this_layer = get_ilayer_from_ID( ID );
  int this_eta = get_ieta_from_ID( ID );
  int this_phi = get_iphi_from_ID( ID );
  
  std::vector<int> adjacent_towers;
  
  // for both IHCal and OHCal, add adjacent layers in the HCal
  if ( this_layer == 0 || this_layer == 1 ) {
    
    for (int delta_layer = 0; delta_layer <= 1; delta_layer++) {
      for (int delta_eta = -1; delta_eta <= 1; delta_eta++) {
	for (int delta_phi = -1; delta_phi <= 1; delta_phi++) {
	  
	  if ( delta_layer == 0 && delta_eta == 0 && delta_phi == 0 ) continue; // this is the same tower
	  
	  int test_eta = this_eta + delta_eta;
	  if ( test_eta < 0 || test_eta >= _HCAL_NETA ) { continue; } // ignore if at the (eta) edge of calorimeter
	  
	  int test_layer = ( this_layer + delta_layer ) % 2; // wrap around in layer
	  int test_phi = ( this_phi + delta_phi + _HCAL_NPHI ) % _HCAL_NPHI; // wrap around in phi (add 64 to avoid -1)
	  
	  // disallow "corner" adjacency (diagonal in eta/phi plane and in different layer) if this option not enabled
	  if ( !_allow_corner_neighbor && delta_layer == 1 && abs( delta_phi ) == 1 && abs( delta_eta ) == 1 ) {
	    if (Verbosity() > 20) std::cout << "RawClusterBuilderTopo::get_adjacent_towers_by_ID : corner growth not allowed " << std::endl;
	    continue;
	  }
	  
	  // add to list of adjacent towers
	  adjacent_towers.push_back( get_ID( test_layer, test_eta, test_phi ) );
	  
	}
      }
    }
  }
  
  // for IHCal only, also add 4x4 group of EMCal towers
  if ( this_layer == 0 && _enable_EMCal ) {
    
    int EMCal_phi_start = get_first_matching_EMCal_phi_from_IHCal( this_phi );
    int EMCal_eta_start = RawClusterBuilderTopo_constants_EMCal_eta_start_given_IHCal[ this_eta ];
    int EMCal_eta_end = RawClusterBuilderTopo_constants_EMCal_eta_end_given_IHCal[ this_eta ];
    
    for (int new_eta = EMCal_eta_start; new_eta <= EMCal_eta_end; new_eta++) {
      for (int delta_phi = 0; delta_phi < 4; delta_phi++) {
	int new_phi = ( EMCal_phi_start + delta_phi + _EMCAL_NPHI ) % _EMCAL_NPHI;
	
	int EMCal_tower = get_ID( 2, new_eta, new_phi );
	if ( Verbosity() > 20 ) std::cout << "RawClusterBuilderTopo::get_adjacent_towers_by_ID : HCal tower with eta / phi = " << this_eta << " / " << this_phi << ", adding EMCal tower with eta / phi = " << new_eta << " / " << new_phi << std::endl;
	adjacent_towers.push_back( EMCal_tower );	
      }
    }
  }
  
  // for EMCal, add adjacent EMCal towers and (usually) one IHCal tower
  if ( this_layer == 2 ) {
    
    // EMCal towers first
    for (int delta_eta = -1; delta_eta <= 1; delta_eta++) {
      for (int delta_phi = -1; delta_phi <= 1; delta_phi++) {
	
	if ( delta_eta == 0 && delta_phi == 0 ) continue; // this is the same tower
	
	int test_eta = this_eta + delta_eta;
	if ( test_eta < 0 || test_eta >= _EMCAL_NETA ) { continue; } // ignore if at the (eta) edge of calorimeter
	
	int test_phi = ( this_phi + delta_phi + _EMCAL_NPHI ) % _EMCAL_NPHI; // wrap around in phi (add 256 to avoid -1)
	
	// add to list of adjacent towers
	adjacent_towers.push_back( get_ID( this_layer, test_eta, test_phi ) );
	
      }
    }
    
    // now add IHCal towers
    if ( _enable_HCal ) {
      int HCal_eta = RawClusterBuilderTopo_constants_IHCal_eta_given_EMCal[ this_eta ];
      int HCal_phi = get_matching_HCal_phi_from_EMCal( this_phi );
      
      if ( HCal_eta >= 0 ) {
	int IHCal_tower = get_ID( 0, HCal_eta, HCal_phi );
	if ( Verbosity() > 20 ) std::cout << "RawClusterBuilderTopo::get_adjacent_towers_by_ID : EMCal tower with eta / phi = " << this_eta << " / " << this_phi << ", adding IHCal tower with eta / phi = " << HCal_eta << " / " << HCal_phi << std::endl;
	adjacent_towers.push_back( IHCal_tower );
      } else {
	if ( Verbosity() > 20 ) std::cout << "RawClusterBuilderTopo::get_adjacent_towers_by_ID : EMCal tower with eta / phi = " << this_eta << " / " << this_phi << ", does not have matching IHCal due to large eta " << std::endl;
      }
    }
    
  }
  
  return adjacent_towers;
  
}

RawClusterBuilderTopo::RawClusterBuilderTopo(const std::string &name)
  : SubsysReco(name)
  , _clusters(nullptr)
{

  // geometry defined at run-time
  _EMCAL_NETA = -1;
  _EMCAL_NPHI = -1;

  _HCAL_NETA = -1;
  _HCAL_NPHI = -1;

  _noise_LAYER[0] = 0.0025;
  _noise_LAYER[1] = 0.006;
  _noise_LAYER[2] = 0.03; // EM

  _sigma_seed = 4.0;
  _sigma_grow = 2.0;
  _sigma_peri = 0.0;

  _allow_corner_neighbor = true;

  _enable_HCal = true;
  _enable_EMCal = true;
  
}

int RawClusterBuilderTopo::InitRun(PHCompositeNode *topNode)
{
  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << PHWHERE << ": " << e.what() << std::endl;
    throw;
  }

  if ( Verbosity() > 0 ) {
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with EMCal enable = " << _enable_EMCal << " and I+OHCal enable = " << _enable_HCal << std::endl;
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with sigma_noise in EMCal / IHCal / OHCal = " << _noise_LAYER[2] << " / " << _noise_LAYER[0] << " / " << _noise_LAYER[1] << std::endl;
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with noise multiples for seeding / growth / perimeter ( S / N / P ) = " << _sigma_seed << " / " << _sigma_grow << " / " << _sigma_peri << std::endl;
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with allow_corner_neighbor = " << _allow_corner_neighbor << " (in HCal)" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterBuilderTopo::process_event(PHCompositeNode *topNode)
{

  RawTowerContainer *towersEM = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  RawTowerContainer *towersIH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");

  if ( !towersEM ) {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWER_CALIB_CEMC does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if ( !towersIH ) {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWER_CALIB_HCALIN does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if ( !towersOH ) {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWER_CALIB_HCALOUT does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  if ( !geomEM ) {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERGEOM_CEMC does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if ( !geomIH ) {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERGEOM_HCALIN does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if ( !geomOH ) {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERGEOM_HCALOUT does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (Verbosity() > 10)
    {
      std::cout << "RawClusterBuilderTopo::process_event: " << towersEM->size() << " TOWER_CALIB_CEMC towers" << std::endl;
      std::cout << "RawClusterBuilderTopo::process_event: " << towersIH->size() << " TOWER_CALIB_HCALIN towers" << std::endl;
      std::cout << "RawClusterBuilderTopo::process_event: " << towersOH->size() << " TOWER_CALIB_HCALOUT towers" << std::endl;

      std::cout << "RawClusterBuilderTopo::process_event: pointer to TOWERGEOM_CEMC: " << geomIH << std::endl;
      std::cout << "RawClusterBuilderTopo::process_event: pointer to TOWERGEOM_HCALIN: " << geomIH << std::endl;
      std::cout << "RawClusterBuilderTopo::process_event: pointer to TOWERGEOM_HCALOUT: " << geomOH << std::endl;
    }


  if ( _EMCAL_NETA < 0 ) {

    // define geometry only once if it has not been yet
    _EMCAL_NETA = geomEM->get_etabins();
    _EMCAL_NPHI = geomEM->get_phibins();

    _EMTOWERMAP_STATUS_ETA_PHI.resize( _EMCAL_NETA, std::vector<int>( _EMCAL_NPHI, -2 ) );
    _EMTOWERMAP_E_ETA_PHI.resize( _EMCAL_NETA, std::vector<float>( _EMCAL_NPHI, 0 ) );
    
  }

  if ( _HCAL_NETA < 0 ) {

    // define geometry only once if it has not been yet
    _HCAL_NETA = geomOH->get_etabins();
    _HCAL_NPHI = geomOH->get_phibins();

    _TOWERMAP_STATUS_LAYER_ETA_PHI.resize( 2, std::vector<std::vector<int> > ( _HCAL_NETA, std::vector<int>( _HCAL_NPHI, -2 ) ) );
    _TOWERMAP_E_LAYER_ETA_PHI.resize( 2, std::vector<std::vector<float> > ( _HCAL_NETA, std::vector<float>( _HCAL_NPHI, 0 ) ) );
    
  }

  // reset maps
  for (int ieta = 0; ieta < _EMCAL_NETA; ieta++) {
    for (int iphi = 0; iphi < _EMCAL_NPHI; iphi++) { 
      
      _EMTOWERMAP_STATUS_ETA_PHI[ ieta ][ iphi ] = -2; // set tower does not exist
      _EMTOWERMAP_E_ETA_PHI[ ieta ][ iphi ] = 0; // set zero energy
      
    }
  }
  for (int ilayer = 0; ilayer < 2; ilayer++) {
    for (int ieta = 0; ieta < _HCAL_NETA; ieta++) {
      for (int iphi = 0; iphi < _HCAL_NPHI; iphi++) { 
	
	_TOWERMAP_STATUS_LAYER_ETA_PHI[ ilayer ][ ieta ][ iphi ] = -2; // set tower does not exist
	_TOWERMAP_E_LAYER_ETA_PHI[ ilayer ][ ieta ][ iphi ] = 0; // set zero energy
	
      }
    }
  }
  
  // setup 
  std::vector< std::pair<int, float> > list_of_seeds;

  // translate towers to our internal representation
  if ( _enable_EMCal ) {
    RawTowerContainer::ConstRange begin_end_EM = towersEM->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
      {
	
	RawTower *tower = rtiter->second;
	RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());
	
	int ieta = geomEM->get_etabin( tower_geom->get_eta() );
	int iphi = geomEM->get_phibin( tower_geom->get_phi() );
	float this_E = tower->get_energy();
	
	_EMTOWERMAP_STATUS_ETA_PHI[ ieta ][ iphi ] = -1; // change status to unknown
	_EMTOWERMAP_E_ETA_PHI[ ieta ][ iphi ] = this_E;
	
	if ( this_E > _sigma_seed * _noise_LAYER[2] ) {
	  int ID = get_ID( 2, ieta, iphi );
	  list_of_seeds.push_back( std::pair<int, float>( ID, this_E ) ); 
	  if (Verbosity() > 10) { 
	    std::cout << "RawClusterBuilderTopo::process_event: adding EMCal tower at ieta / iphi = " << ieta << " / " << iphi << " with E = " << this_E << std::endl;
	    std::cout << " --> ID = " << ID << " , check ilayer / ieta / iphi = " << get_ilayer_from_ID( ID ) << " / " << get_ieta_from_ID( ID ) << " / " << get_iphi_from_ID( ID ) << std::endl;
	  };
	}
	
      }
  }

  // translate towers to our internal representation
  if ( _enable_HCal ) {
    RawTowerContainer::ConstRange begin_end_IH = towersIH->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_IH.first; rtiter != begin_end_IH.second; ++rtiter)
      {
	
	RawTower *tower = rtiter->second;
	RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
	
	int ieta = geomIH->get_etabin( tower_geom->get_eta() );
	int iphi = geomIH->get_phibin( tower_geom->get_phi() );
	float this_E = tower->get_energy();
	
	_TOWERMAP_STATUS_LAYER_ETA_PHI[ 0 ][ ieta ][ iphi ] = -1; // change status to unknown
	_TOWERMAP_E_LAYER_ETA_PHI[ 0 ][ ieta ][ iphi ] = this_E;
	
	if ( this_E > _sigma_seed * _noise_LAYER[0] ) {
	  int ID = get_ID( 0, ieta, iphi );
	  list_of_seeds.push_back( std::pair<int, float>( ID, this_E ) ); 
	  if (Verbosity() > 10) { 
	    std::cout << "RawClusterBuilderTopo::process_event: adding IHCal tower at ieta / iphi = " << ieta << " / " << iphi << " with E = " << this_E << std::endl;
	    std::cout << " --> ID = " << ID << " , check ilayer / ieta / iphi = " << get_ilayer_from_ID( ID ) << " / " << get_ieta_from_ID( ID ) << " / " << get_iphi_from_ID( ID ) << std::endl;
	  };
	}
	
      }
    
    RawTowerContainer::ConstRange begin_end_OH = towersOH->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_OH.first; rtiter != begin_end_OH.second; ++rtiter)
      {
	
	RawTower *tower = rtiter->second;
	RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());
	
	int ieta = geomOH->get_etabin( tower_geom->get_eta() );
	int iphi = geomOH->get_phibin( tower_geom->get_phi() );
	float this_E = tower->get_energy();
	
	_TOWERMAP_STATUS_LAYER_ETA_PHI[ 1 ][ ieta ][ iphi ] = -1; // change status to unknown
	_TOWERMAP_E_LAYER_ETA_PHI[ 1 ][ ieta ][ iphi ] = this_E;
	
	if ( this_E > _sigma_seed * _noise_LAYER[1] ) {
	  int ID = get_ID( 1, ieta, iphi );
	  list_of_seeds.push_back( std::pair<int, float>( ID, this_E ) ); 
	  if (Verbosity() > 10) { 
	    std::cout << "RawClusterBuilderTopo::process_event: adding OHCal tower at ieta / iphi = " << ieta << " / " << iphi << " with E = " << this_E << std::endl;
	    std::cout << " --> ID = " << ID << " , check ilayer / ieta / iphi = " << get_ilayer_from_ID( ID ) << " / " << get_ieta_from_ID( ID ) << " / " << get_iphi_from_ID( ID ) << std::endl;
	  };
	}
	
      }
  }
  
  if (Verbosity() > 10) {
    for (unsigned int n = 0; n < list_of_seeds.size(); n++) {
      std::cout << "RawClusterBuilderTopo::process_event: unsorted seed element n = " << n << " , ID / E = " << list_of_seeds.at(n).first << " / " << list_of_seeds.at(n).second << std::endl;
    }
  }

  std::sort( list_of_seeds.begin(), list_of_seeds.end(), sort_by_pair_second );

  if (Verbosity() > 10) {
    for (unsigned int n = 0; n < list_of_seeds.size(); n++) {
      std::cout << "RawClusterBuilderTopo::process_event: sorted seed element n = " << n << " , ID / E = " << list_of_seeds.at(n).first << " / " << list_of_seeds.at(n).second << std::endl;
    }
  }

  if (Verbosity() > 0)
    std::cout << "RawClusterBuilderTopo::process_event: initialized with " << list_of_seeds.size() << " seeds with E > 4*sigma " << std::endl;


  int cluster_index = 0; // begin counting clusters

  while ( list_of_seeds.size() > 0 ) {
   
    int seed_ID = list_of_seeds.at( 0 ).first;
    list_of_seeds.erase( list_of_seeds.begin() );

    if (Verbosity() > 5) {
      std::cout << " RawClusterBuilderTopo::process_event: in seeded loop, current seed has ID = " << seed_ID << " , length of remaining seed vector = " << list_of_seeds.size() << std::endl;
    }
    
    // if this seed was already claimed by some other seed during its growth, remove it and do nothing
    int seed_status = get_status_from_ID( seed_ID );
    if ( seed_status > -1 ) {
      if (Verbosity() > 10)
	std::cout << " --> already owned by cluster # " << seed_status << std::endl;
      continue; // go onto the next iteration of the loop
    }

    // this seed tower now owned by new cluster
    set_status_by_ID( seed_ID, cluster_index );

    std::vector<int> cluster_tower_ID;
    cluster_tower_ID.push_back( seed_ID );

    std::vector<int> grow_tower_ID;
    grow_tower_ID.push_back( seed_ID );

    // iteratively process growth towers, adding > 2 * sigma neighbors to the list for further checking

    if (Verbosity() > 5)
      std::cout << " RawClusterBuilderTopo::process_event: Entering Growth stage for cluster " << cluster_index << std::endl;
    
    while ( grow_tower_ID.size() > 0 ) {

      int grow_ID = grow_tower_ID.at( 0 );
      grow_tower_ID.erase( grow_tower_ID.begin() );
      
      if (Verbosity() > 5)
	std::cout << " --> cluster " << cluster_index << ", growth stage, examining neighbors of ID " << grow_ID << ", " << grow_tower_ID.size() << " grow towers left" << std::endl;

      std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID( grow_ID );

      for ( unsigned int adj = 0; adj < adjacent_tower_IDs.size(); adj++) {

	int this_adjacent_tower_ID = adjacent_tower_IDs.at( adj );
	
	if (Verbosity() > 10) std::cout << " --> --> --> checking possible adjacent tower with ID " << this_adjacent_tower_ID << " : ";

	int test_layer = get_ilayer_from_ID( this_adjacent_tower_ID );
	
	// if tower does not exist, continue
	if ( get_status_from_ID( this_adjacent_tower_ID ) == -2 ) {
	  if (Verbosity() > 10) std::cout << "does not exist " << std::endl;
	  continue;
	}

	// if tower is owned by THIS cluster already, continue
	if ( get_status_from_ID( this_adjacent_tower_ID ) == cluster_index ) {
	  if (Verbosity() > 10) std::cout << "already owned by this cluster index " << cluster_index << std::endl;
	  continue;
	}
	    
	// if tower has < 2*sigma energy, continue
	if ( get_E_from_ID( this_adjacent_tower_ID ) < _sigma_grow * _noise_LAYER[ test_layer ] ) {
	  if (Verbosity() > 10) std::cout << "E = " << get_E_from_ID( this_adjacent_tower_ID ) << " under 2*sigma threshold " << std::endl;
	  continue;
	}

	// if tower is owned by somebody else, continue (although should this really happen?)
	if ( get_status_from_ID( this_adjacent_tower_ID ) > -1  ) {
	  if (Verbosity() > 10) std::cout << "ERROR! in growth stage, encountered >2sigma tower which is already owned?!" << std::endl;
	  continue;
	}
	  
	// tower good to be added to cluster and to list of grow towers
	grow_tower_ID.push_back( this_adjacent_tower_ID );
	cluster_tower_ID.push_back( this_adjacent_tower_ID );
	set_status_by_ID( this_adjacent_tower_ID, cluster_index );
	if (Verbosity() > 10) std::cout << "add this tower ( ID " <<  this_adjacent_tower_ID << " ) to grow list " << std::endl;
	
      }

      if (Verbosity() > 5) std::cout << " --> after examining neighbors, grow list is now " <<  grow_tower_ID.size() << ", # of towers in cluster = " << cluster_tower_ID.size() << std::endl;

    }

    // done growing cluster, now add on perimeter towers with E > 0 * sigma
    if (Verbosity() > 5)
      std::cout << " RawClusterBuilderTopo::process_event: Entering Perimeter stage for cluster " << cluster_index << std::endl;

    // we'll be adding on to the cluster list, so get the # of core towers first 
    int n_core_towers = cluster_tower_ID.size();

    for ( int ic = 0; ic < n_core_towers; ic++) {
      
      int core_ID = cluster_tower_ID.at( ic );
      
      if (Verbosity() > 5)
	std::cout << " --> cluster " << cluster_index << ", perimeter stage, examining neighbors of ID " << core_ID << ", core cluster # " << ic << " of " << n_core_towers << " total " << std::endl;

      std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID( core_ID );

      for ( unsigned int adj = 0; adj < adjacent_tower_IDs.size(); adj++) {

	int this_adjacent_tower_ID = adjacent_tower_IDs.at( adj );

	if (Verbosity() > 10) std::cout << " --> --> --> checking possible adjacent tower with ID " << this_adjacent_tower_ID << " : ";

	int test_layer = get_ilayer_from_ID( this_adjacent_tower_ID );
	
	// if tower does not exist, continue
	if ( get_status_from_ID( this_adjacent_tower_ID ) == -2 ) {
	  if (Verbosity() > 10) std::cout << "does not exist " << std::endl;
	  continue;
	}
	
	// if tower is owned by somebody else (including current cluster), continue. ( allowed during perimeter fixing state )
	if ( get_status_from_ID( this_adjacent_tower_ID ) > -1  ) {
	  if (Verbosity() > 10) std::cout << "already owned by other cluster index " <<  get_status_from_ID( this_adjacent_tower_ID ) << std::endl;
	  continue;
	}
	
	// if tower has < 0*sigma energy, continue
	if ( get_E_from_ID( this_adjacent_tower_ID ) < _sigma_peri * _noise_LAYER[ test_layer ] ) {
	  if (Verbosity() > 10) std::cout << "E = " << get_E_from_ID( this_adjacent_tower_ID ) << " under 0*sigma threshold " << std::endl;
	  continue;
	}
	
	// perimeter tower good to be added to cluster
	cluster_tower_ID.push_back( this_adjacent_tower_ID );
	set_status_by_ID( this_adjacent_tower_ID, cluster_index );
	if (Verbosity() > 10) std::cout << "add this tower ( ID " <<  this_adjacent_tower_ID << " ) to cluster " << std::endl;
	
      }
      
      if (Verbosity() > 5) std::cout << " --> after examining perimeter neighbors, # of towers in cluster is now = " << cluster_tower_ID.size() << std::endl;
    }

    // increment cluster index for next one
    cluster_index++;
    
  }

  if (Verbosity() > 0) std::cout << "RawClusterBuilderTopo::process_event: " << cluster_index << " topo-clusters reconstructed, translating to cluster form:" << std::endl;

  // defnite temporary objects here
  std::vector<RawCluster*> clusters;
  std::vector<float> cluster_E;
  std::vector<float> cluster_x;
  std::vector<float> cluster_y;
  std::vector<float> cluster_z;

  for (int cl = 0; cl < cluster_index; cl++) {
    clusters.push_back( new RawClusterv1() );
    cluster_E.push_back( 0 );
    cluster_x.push_back( 0 );
    cluster_y.push_back( 0 );
    cluster_z.push_back( 0 );
  }

  // fill these using original tower information & internal grid 

  {
    RawTowerContainer::ConstRange begin_end_EM = towersEM->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
      {
	
	RawTower *tower = rtiter->second;
	RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());
	
	int ieta = geomEM->get_etabin( tower_geom->get_eta() );
	int iphi = geomEM->get_phibin( tower_geom->get_phi() );
	float this_E = tower->get_energy();
	
	int status = _EMTOWERMAP_STATUS_ETA_PHI[ ieta ][ iphi ];
	if ( status < 0 ) continue;
	
	clusters[ status ]->addTower( tower->get_id(), this_E );
	cluster_E[ status ] = cluster_E[ status ] + this_E ;
	cluster_x[ status ] = cluster_x[ status ] + this_E * tower_geom->get_center_x();
	cluster_y[ status ] = cluster_y[ status ] + this_E * tower_geom->get_center_y();
	cluster_z[ status ] = cluster_z[ status ] + this_E * tower_geom->get_center_z();
	
      }
  }

  {
    RawTowerContainer::ConstRange begin_end_IH = towersIH->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_IH.first; rtiter != begin_end_IH.second; ++rtiter)
      {
	
	RawTower *tower = rtiter->second;
	RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
	
	int ieta = geomIH->get_etabin( tower_geom->get_eta() );
	int iphi = geomIH->get_phibin( tower_geom->get_phi() );
	float this_E = tower->get_energy();
	
	int status = _TOWERMAP_STATUS_LAYER_ETA_PHI[ 0 ][ ieta ][ iphi ];
	if ( status < 0 ) continue;
	
	clusters[ status ]->addTower( tower->get_id(), this_E );
	cluster_E[ status ] = cluster_E[ status ] + this_E ;
	cluster_x[ status ] = cluster_x[ status ] + this_E * tower_geom->get_center_x();
	cluster_y[ status ] = cluster_y[ status ] + this_E * tower_geom->get_center_y();
	cluster_z[ status ] = cluster_z[ status ] + this_E * tower_geom->get_center_z();
	
      }
  }

  {
    RawTowerContainer::ConstRange begin_end_OH = towersOH->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_OH.first; rtiter != begin_end_OH.second; ++rtiter)
      {
	
	RawTower *tower = rtiter->second;
	RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());
	
	int ieta = geomOH->get_etabin( tower_geom->get_eta() );
	int iphi = geomOH->get_phibin( tower_geom->get_phi() );
	float this_E = tower->get_energy();
	
	int status = _TOWERMAP_STATUS_LAYER_ETA_PHI[ 1 ][ ieta ][ iphi ];
	if ( status < 0 ) continue;
	
	clusters[ status ]->addTower( tower->get_id(), this_E );
	cluster_E[ status ] = cluster_E[ status ] + this_E ;
	cluster_x[ status ] = cluster_x[ status ] + this_E * tower_geom->get_center_x();
	cluster_y[ status ] = cluster_y[ status ] + this_E * tower_geom->get_center_y();
	cluster_z[ status ] = cluster_z[ status ] + this_E * tower_geom->get_center_z();
	
      }
  }

  // iterate through and add to official container

  for (int cl = 0; cl < cluster_index; cl++) {
    
    clusters[ cl ]->set_energy( cluster_E[ cl ] );
    
    float mean_x = cluster_x[ cl ] /  cluster_E[ cl ];
    float mean_y = cluster_y[ cl ] /  cluster_E[ cl ];
    float mean_z = cluster_z[ cl ] /  cluster_E[ cl ];

    clusters[ cl ]->set_r( sqrt( mean_y * mean_y + mean_x * mean_x) );
    clusters[ cl ]->set_phi( atan2( mean_y, mean_x) );
    clusters[ cl ]->set_z( mean_z );
    
    _clusters->AddCluster( clusters[ cl ] );

  }
  
  if ( Verbosity() > 0 ) {
    RawClusterContainer::ConstRange begin_end = _clusters->getClusters();
    for (RawClusterContainer::ConstIterator hiter = begin_end.first; hiter != begin_end.second; ++hiter)
      {
	std::cout << "--> ";
	hiter->second->identify();
	std::cout << std::endl;
      }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterBuilderTopo::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderTopo::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in RawClusterBuilderTopo::CreateNodes");
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TOPOCLUSTER"));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode("TOPOCLUSTER");
    dstNode->addNode(DetNode);
  }

  _clusters = new RawClusterContainer();
  ClusterNodeName = "TOPOCLUSTER_HCAL";
  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, ClusterNodeName.c_str(), "PHObject");
  DetNode->addNode(clusterNode);
}
