#include "DetermineTowerBackground.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>

// sPHENIX includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include "TowerBackground_v1.h"

#include <TMath.h>

// standard includes
#include <iomanip>
#include <iostream>
#include <vector>

DetermineTowerBackground::DetermineTowerBackground(const std::string &name)
  : SubsysReco(name)
{

  _do_flow = false;

  _seed_jet_D = 3.0;
  _seed_jet_pt = 7.0;

  _v2 = 0;
  _Psi2 = 0;

  _UE.resize(3, std::vector<float>(1, 0) );

  // initiate sizes as -1 to tell module they should be set when it
  // sees the HCal geometry for the first time

  _HCAL_NETA = -1;
  _HCAL_NPHI = -1;

  _backgroundName = "TestTowerBackground";
  _seed_type = 0;
}

DetermineTowerBackground::~DetermineTowerBackground()
{
}


void DetermineTowerBackground::SetBackgroundOutputName( std::string name ) {
  
  _backgroundName = name;
  
  return;
}

void DetermineTowerBackground::SetSeedType( int seed_type ) {

  _seed_type = seed_type;

  return;

}

int DetermineTowerBackground::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    std::cout << "DetermineTowerBackground::Init: initialized" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int DetermineTowerBackground::InitRun(PHCompositeNode *topNode)
{
  return CreateNode(topNode);
}

int DetermineTowerBackground::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0) {
    std::cout << "DetermineTowerBackground::process_event: entering with do_flow = " << _do_flow << ", seed type = " << _seed_type << ", ";
    if ( _seed_type == 0 ) std::cout << " D = " << _seed_jet_D << std::endl;
    else if ( _seed_type == 1 ) std::cout << " pT = " << _seed_jet_pt << std::endl;
    else std::cout << " UNDEFINED seed behavior! " << std::endl;
  }

  // clear seed eta/phi positions
  _seed_eta.resize(0);
  _seed_phi.resize(0);

  // pull out the tower containers and geometry objects at the start
  RawTowerContainer *towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
  RawTowerContainer *towersIH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");
  if (Verbosity() > 0) {
    std::cout << "DetermineTowerBackground::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC_RETOWER towers" << std::endl;
    std::cout << "DetermineTowerBackground::process_event: " << towersIH3->size() << " TOWER_CALIB_HCALIN towers" << std::endl;
    std::cout << "DetermineTowerBackground::process_event: " << towersOH3->size() << " TOWER_CALIB_HCALOUT towers" << std::endl;
  }

  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  // seed type 0 is D > 3 R=0.2 jets run on retowerized CEMC
  if (_seed_type == 0) {
    JetMap* reco2_jets = findNode::getClass<JetMap>(topNode,"AntiKt_Tower_HIRecoSeedsRaw_r02");

    if (Verbosity() > 1)
      std::cout << "DetermineTowerBackground::process_event: examining possible seeds (1st iteration) ... " << std::endl;

    for (JetMap::Iter iter = reco2_jets->begin(); iter != reco2_jets->end(); ++iter) {

      Jet* this_jet = iter->second;

      float this_pt = this_jet->get_pt();
      float this_phi = this_jet->get_phi();
      float this_eta = this_jet->get_eta();

      if (this_jet->get_pt() < 5) continue;

      if (Verbosity() > 2)
	std::cout << "DetermineTowerBackground::process_event: possible seed jet with pt / eta / phi = " << this_pt << " / " << this_eta << " / " << this_phi << ", examining constituents..." << std::endl;

      std::map< int, double > constituent_ETsum;

      for (Jet::ConstIter comp = this_jet->begin_comp(); comp !=  this_jet->end_comp(); ++comp) {

	int comp_ieta = -1;
	int comp_iphi = -1;
	float comp_ET = 0;

	RawTower *tower;
	RawTowerGeom *tower_geom;
	
        if ( (*comp).first == 5 ) {
          tower = towersIH3->getTower( (*comp).second );
          tower_geom = geomIH->get_tower_geometry(tower->get_key());

	  comp_ieta = geomIH->get_etabin( tower_geom->get_eta() );
	  comp_iphi =  geomIH->get_phibin( tower_geom->get_phi() );
	  comp_ET = tower->get_energy() / cosh( tower_geom->get_eta() );
        }
        else if ( (*comp).first == 7 ) {
          tower = towersOH3->getTower( (*comp).second );
          tower_geom = geomOH->get_tower_geometry(tower->get_key());

	  comp_ieta = geomIH->get_etabin( tower_geom->get_eta() );
	  comp_iphi =  geomIH->get_phibin( tower_geom->get_phi() );
	  comp_ET = tower->get_energy() / cosh( tower_geom->get_eta() );
        }
        else if ( (*comp).first == 13 ) {
          tower = towersEM3->getTower( (*comp).second );
          tower_geom = geomIH->get_tower_geometry(tower->get_key());

	  comp_ieta = geomIH->get_etabin( tower_geom->get_eta() );
	  comp_iphi =  geomIH->get_phibin( tower_geom->get_phi() );
	  comp_ET = tower->get_energy() / cosh( tower_geom->get_eta() );
	}

	int comp_ikey = 1000 * comp_ieta + comp_iphi;

	if (Verbosity() > 4)
	  std::cout << "DetermineTowerBackground::process_event: --> --> constituent in layer " << (*comp).first << " at ieta / iphi = " << comp_ieta << " / " << comp_iphi << ", filling map with key = " << comp_ikey << " and ET = " << comp_ET << std::endl;

	constituent_ETsum[ comp_ikey ] += comp_ET;

	if (Verbosity() > 4)
	  std::cout << "DetermineTowerBackground::process_event: --> --> ET sum map at key = " << comp_ikey << " now has ET = " << constituent_ETsum[ comp_ikey ] << std::endl;
	
      }

      // now iterate over constituent_ET sums to find maximum and mean
      float constituent_max_ET = 0;
      float constituent_sum_ET = 0;
      int nconstituents = 0;
      
      if (Verbosity() > 4)
	std::cout << "DetermineTowerBackground::process_event: --> now iterating over map..." << std::endl;
      for (std::map<int,double>::iterator map_iter = constituent_ETsum.begin(); map_iter != constituent_ETsum.end(); ++map_iter) {
	if (Verbosity() > 4)
	  std::cout << "DetermineTowerBackground::process_event: --> --> map has key # " << map_iter->first << " and ET = " << map_iter->second << std::endl;
	nconstituents++;
	constituent_sum_ET +=  map_iter->second;
	if ( map_iter->second > constituent_max_ET ) constituent_max_ET = map_iter->second;
      }

      float mean_constituent_ET = constituent_sum_ET / nconstituents;
      float seed_D = constituent_max_ET / mean_constituent_ET;
      
      if (Verbosity() > 3)
	std::cout << "DetermineTowerBackground::process_event: --> jet has < ET > = " << constituent_sum_ET << " / " << nconstituents << " = " << mean_constituent_ET << ", max-ET = " << constituent_max_ET << ", and D = " << seed_D << std::endl;
      
      if ( seed_D > _seed_jet_D ) {
	_seed_eta.push_back( this_eta );
	_seed_phi.push_back( this_phi );
	
	if (Verbosity() > 1)
	  std::cout << "DetermineTowerBackground::process_event: --> adding seed at eta / phi = " << this_eta << " / " << this_phi << " ( R=0.2 jet with pt = " << this_pt << ", D = " << seed_D << " ) " << std::endl;
      } else {
	if (Verbosity() > 3)
	  std::cout << "DetermineTowerBackground::process_event: --> discarding potential seed at eta / phi = " << this_eta << " / " << this_phi << " ( R=0.2 jet with pt = " << this_pt << ", D = " << seed_D << " ) " << std::endl;
      }
      
    }

  }

  // seed type 1 is the set of those jets above which, when their
  // kinematics are updated for the first background subtraction, have
  // pT > 20 GeV
  if (_seed_type == 1) {
    JetMap* reco2_jets = findNode::getClass<JetMap>(topNode,"AntiKt_Tower_HIRecoSeedsSub_r02");

    if (Verbosity() > 1)
      std::cout << "DetermineTowerBackground::process_event: examining possible seeds (2nd iteration) ... " << std::endl;
    
    for (JetMap::Iter iter = reco2_jets->begin(); iter != reco2_jets->end(); ++iter) {
      Jet* this_jet = iter->second;
      
      float this_pt = this_jet->get_pt();
      float this_phi = this_jet->get_phi();
      float this_eta = this_jet->get_eta();

      if (this_jet->get_pt() < _seed_jet_pt ) continue;

      _seed_eta.push_back( this_eta );
      _seed_phi.push_back( this_phi );

      if (Verbosity() > 1)
	std::cout << "DetermineTowerBackground::process_event: --> adding seed at eta / phi = " << this_eta << " / " << this_phi << " ( R=0.2 jet with pt = " << this_pt << " ) " << std::endl;
    }
    
  }


  // get the binning from the geometry (different for 1D vs 2D...)
  if (_HCAL_NETA < 0) {
    
    _HCAL_NETA = geomIH->get_etabins();
    _HCAL_NPHI = geomIH->get_phibins();
    
    // resize UE density and energy vectors 
    _UE[0].resize( _HCAL_NETA, 0 );
    _UE[1].resize( _HCAL_NETA, 0 );
    _UE[2].resize( _HCAL_NETA, 0 );

    _EMCAL_E.resize(_HCAL_NETA, std::vector<float>(_HCAL_NPHI, 0));
    _IHCAL_E.resize(_HCAL_NETA, std::vector<float>(_HCAL_NPHI, 0));
    _OHCAL_E.resize(_HCAL_NETA, std::vector<float>(_HCAL_NPHI, 0));

    // for flow determination, build up a 1-D phi distribution of
    // energies from all layers summed together, populated only from eta
    // strips which do not have any excluded phi towers
    _FULLCALOFLOW_PHI_E.resize( _HCAL_NPHI, 0 );
    _FULLCALOFLOW_PHI_VAL.resize( _HCAL_NPHI, 0 );

    if (Verbosity() > 0) {
      std::cout << "DetermineTowerBackground::process_event: setting number of towers in eta / phi: " << _HCAL_NETA << " / " << _HCAL_NPHI << std::endl;
    }

  }

  // reset all maps map
  for (int ieta = 0; ieta < _HCAL_NETA; ieta++) {
    for (int iphi = 0; iphi < _HCAL_NPHI; iphi++) {
      _EMCAL_E[ieta][iphi] = 0;
      _IHCAL_E[ieta][iphi] = 0;
      _OHCAL_E[ieta][iphi] = 0;
    }
  }

  for (int iphi = 0; iphi < _HCAL_NPHI; iphi++) {
    _FULLCALOFLOW_PHI_E[iphi] = 0;
    _FULLCALOFLOW_PHI_VAL[iphi] = 0;
  }
  
  // iterate over EMCal towers
  RawTowerContainer::ConstRange begin_end = towersEM3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter) {
    RawTower *tower = rtiter->second;

    RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
    
    float this_eta = tower_geom->get_eta();
    float this_phi = tower_geom->get_phi();
    int this_etabin = geomIH->get_etabin(this_eta);
    int this_phibin = geomIH->get_phibin(this_phi);
    float this_E = tower->get_energy();

    _EMCAL_E[ this_etabin ][ this_phibin ] += this_E;
    
    if (Verbosity() > 2 && tower->get_energy() > 1)
      {
	std::cout << "DetermineTowerBackground::process_event: EMCal tower eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << this_etabin << " ) / " << this_phi << " ( " << this_phibin << " ) / " << this_E << std::endl;
      }
  }

  // iterate over IHCal towers
  RawTowerContainer::ConstRange begin_end_IH = towersIH3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_IH.first; rtiter != begin_end_IH.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());

    float this_eta = tower_geom->get_eta();
    float this_phi = tower_geom->get_phi();
    int this_etabin = geomIH->get_etabin( this_eta );
    int this_phibin = geomIH->get_phibin( this_phi );
    float this_E = tower->get_energy();

    _IHCAL_E[ this_etabin ][ this_phibin ] += this_E;

    if (Verbosity() > 2 && tower->get_energy() > 1)
    {
      std::cout << "DetermineTowerBackground::process_event: IHCal tower at eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << this_etabin << " ) / " << this_phi << " ( " << this_phibin << " ) / " << this_E << std::endl;
    }
  }

  // iterate over OHCal towers
  RawTowerContainer::ConstRange begin_end_OH = towersOH3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_OH.first; rtiter != begin_end_OH.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());

    float this_eta = tower_geom->get_eta();
    float this_phi = tower_geom->get_phi();
    int this_etabin = geomOH->get_etabin( this_eta );
    int this_phibin = geomOH->get_phibin( this_phi );
    float this_E = tower->get_energy();

    _OHCAL_E[ this_etabin ][ this_phibin ] += this_E;

    if (Verbosity() > 2 && tower->get_energy() > 1)
    {
      std::cout << "DetermineTowerBackground::process_event: OHCal tower at eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << this_etabin << " ) / " << this_phi << " ( " << this_phibin << " ) / " << this_E << std::endl;
    }
  }

  // first, calculate flow: Psi2 & v2, if enabled

  _Psi2 = 0;
  _v2 = 0;

  if ( _do_flow ) {
    // check for the case when every tower is excluded
    int nStripsAvailableForFlow = 0;
    int nStripsUnavailableForFlow = 0;
    
    for (int layer = 0; layer < 3; layer++) {
      
      int local_max_eta = _HCAL_NETA;
      int local_max_phi = _HCAL_NPHI;
      
      for (int eta = 0; eta < local_max_eta; eta++) {
	
	bool isAnyTowerExcluded = false;
	
	for (int phi = 0; phi < local_max_phi; phi++) {
	  float this_eta = geomIH->get_etacenter( eta );
	  float this_phi = geomIH->get_phicenter( phi );
	  
	  bool isExcluded = false;
	  
	  for (unsigned int iseed = 0; iseed < _seed_eta.size(); iseed++) {
	    
	    float deta = this_eta - _seed_eta[ iseed ];
	    float dphi = this_phi - _seed_phi[ iseed ];
	    if (dphi >  3.14159) dphi -= 2*3.14159;
	    if (dphi < -3.14159) dphi += 2*3.14159;
	    float dR = sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
	    if (dR < 0.4) {
	      isExcluded = true;
	      if (Verbosity() > 10) std::cout << " setting excluded mark from seed at eta / phi = " << _seed_eta[ iseed ] << " / " << _seed_phi[ iseed ] << std::endl;
	    }
	  }
	  
	  // if even a single tower in this eta strip is excluded, we
	  // can't use the strip for flow determination
	  if (isExcluded)
	    isAnyTowerExcluded = true;
	} // close phi loop
	
	// if this eta strip can be used for flow determination, fill it now
	if (!isAnyTowerExcluded) {
	  
	  if ( Verbosity() > 4) 
	    std::cout << " strip at layer " << layer << ", eta " << eta << " has no excluded towers and can be used for flow determination " << std::endl;
	  nStripsAvailableForFlow++;
	  
	  for (int phi = 0; phi < local_max_phi; phi++) {
	    float this_phi = geomIH->get_phicenter( phi );
	    
	    if ( layer == 0 ) _FULLCALOFLOW_PHI_E[ phi ] += _EMCAL_E[ eta ][ phi ];
	    if ( layer == 1 ) _FULLCALOFLOW_PHI_E[ phi ] += _IHCAL_E[ eta ][ phi ];
	    if ( layer == 2 ) _FULLCALOFLOW_PHI_E[ phi ] += _OHCAL_E[ eta ][ phi ];
	    
	    _FULLCALOFLOW_PHI_VAL[ phi ] = this_phi; // should really set this globally only one time
	    
	  }
	  
	} else {
	  if ( Verbosity() > 4) 
	    std::cout << " strip at layer " << layer << ", eta " << eta << " DOES have excluded towers and CANNOT be used for flow determination " << std::endl;
	  nStripsUnavailableForFlow++;
	}
	
      } // close eta loop
    } // close layer loop
    
    // flow determination
    
    float Q_x = 0;
    float Q_y = 0;
    float E = 0;
    
    float sum_cos2dphi = 0;
    
    if (Verbosity() > 0 )
      std::cout << "DetermineTowerBackground::process_event: # of strips (summed over layers) available / unavailable for flow determination: " << nStripsAvailableForFlow << " / " << nStripsUnavailableForFlow << std::endl;
    
    if ( nStripsAvailableForFlow > 0 ) { 
      for (int iphi = 0; iphi < _HCAL_NPHI; iphi++) {
	
	E += _FULLCALOFLOW_PHI_E[ iphi ];
	Q_x += _FULLCALOFLOW_PHI_E[ iphi ] * cos( 2 * _FULLCALOFLOW_PHI_VAL[ iphi ] );
	Q_y += _FULLCALOFLOW_PHI_E[ iphi ] * sin( 2 * _FULLCALOFLOW_PHI_VAL[ iphi ] );
	
      }
      
      _Psi2 = TMath::ATan2( Q_y, Q_x ) / 2.0 ;
      
      for (int iphi = 0; iphi < _HCAL_NPHI; iphi++) {
	
	sum_cos2dphi += _FULLCALOFLOW_PHI_E[ iphi ] * cos( 2 * (  _FULLCALOFLOW_PHI_VAL[ iphi ] - _Psi2 ) );
	
      }
      
      _v2 = sum_cos2dphi / E;
      
    } else {
      _Psi2 = 0;
      _v2 = 0;
      if (Verbosity() > 0 )
	std::cout << "DetermineTowerBackground::process_event: no full strips available for flow modulation, setting v2 and Psi = 0" << std::endl;
    }

    if (Verbosity() > 0 ) {
      
      std::cout << "DetermineTowerBackground::process_event: unnormalized Q vector (Qx, Qy) = ( " << Q_x << ", " << Q_y << " ) with Sum E_i = " << E << std::endl;
      std::cout << "DetermineTowerBackground::process_event: Psi2 = " << _Psi2 << " ( " << _Psi2 / 3.14159 << " * pi ) , v2 = " << _v2 << std::endl;
    }
  } // if do flow
  else {
    if (Verbosity() > 0 ) {
      std::cout << "DetermineTowerBackground::process_event: flow not enabled, setting Psi2 = " << _Psi2 << " ( " << _Psi2 / 3.14159 << " * pi ) , v2 = " << _v2 << std::endl;
    }
  }

  // now calculate energy densities...

  // starting with the EMCal first...
  for (int layer = 0; layer < 3; layer++) {
    
    int local_max_eta = _HCAL_NETA;
    int local_max_phi = _HCAL_NPHI;

    for (int eta = 0; eta < local_max_eta; eta++) {

      float total_E = 0;
      int total_tower = 0;

      for (int phi = 0; phi < local_max_phi; phi++) {
	float this_eta = geomIH->get_etacenter( eta );
	float this_phi = geomIH->get_phicenter( phi );
	
	bool isExcluded = false;

	for (unsigned int iseed = 0; iseed < _seed_eta.size(); iseed++) {
	  
	  float deta = this_eta - _seed_eta[ iseed ];
	  float dphi = this_phi - _seed_phi[ iseed ];
	  if (dphi >  3.14159) dphi -= 2*3.14159;
	  if (dphi < -3.14159) dphi += 2*3.14159;
	  float dR = sqrt( pow( deta, 2 ) + pow( dphi, 2 ) );
	  if (dR < 0.4) {
	    isExcluded = true;
	    if (Verbosity() > 10) std::cout << " setting excluded mark from seed at eta / phi = " << _seed_eta[ iseed ] << " / " << _seed_phi[ iseed ] << std::endl;
	  }
	}

	if (!isExcluded) {
	  if ( layer == 0 ) total_E += _EMCAL_E[ eta ][ phi ] / ( 1 + 2 * _v2 * cos( 2 * ( this_phi - _Psi2 ) ) );
	  if ( layer == 1 ) total_E += _IHCAL_E[ eta ][ phi ] / ( 1 + 2 * _v2 * cos( 2 * ( this_phi - _Psi2 ) ) );
	  if ( layer == 2 ) total_E += _OHCAL_E[ eta ][ phi ] / ( 1 + 2 * _v2 * cos( 2 * ( this_phi - _Psi2 ) ) );
	  total_tower++;
	} else {
	  if (Verbosity() > 10) std::cout << " tower at eta / phi = " << this_eta << " / " << this_phi << " with E = " << total_E << " excluded due to seed " << std::endl;
	}

      }

      std::pair< float, float > etabounds = geomIH->get_etabounds( eta );
      std::pair< float, float > phibounds = geomIH->get_phibounds( 0 );

      float deta = etabounds.second - etabounds.first;
      float dphi = phibounds.second - phibounds.first;
      float total_area = total_tower * deta * dphi;
      _UE[ layer ].at( eta ) = total_E / total_tower;
      
      if (Verbosity() > 3 ) {
	std::cout << "DetermineTowerBackground::process_event: at layer / eta index ( eta range ) = " << layer << " / " << eta << " ( " << etabounds.first << " - " << etabounds.second << " ) , total E / total Ntower / total area = " << total_E << " / " << total_tower << " / " << total_area << " , UE per tower = " << total_E / total_tower << std::endl;
      }
      
    }
  }

  if (Verbosity() > 0 ) {

    for (int layer = 0; layer < 3; layer++) {
      std::cout << "DetermineTowerBackground::process_event: summary UE in layer " << layer << " : ";
      for (int eta = 0; eta < _HCAL_NETA; eta++) std::cout << _UE[ layer ].at( eta ) << " , ";
      std::cout << std::endl;
    }

  }

  // 

  FillNode(topNode);

  if (Verbosity() > 0) std::cout << "DetermineTowerBackground::process_event: exiting" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int DetermineTowerBackground::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int DetermineTowerBackground::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the jet background stuff under a sub-node directory
  PHCompositeNode *bkgNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "JETBACKGROUND"));
  if (!bkgNode)
  {
    bkgNode = new PHCompositeNode("JETBACKGROUND");
    dstNode->addNode(bkgNode);
  }

  // create the TowerBackground node...
  TowerBackground *towerbackground = findNode::getClass<TowerBackground>(topNode, _backgroundName );
  if (!towerbackground)
  {
    towerbackground = new TowerBackground_v1();
    PHIODataNode<PHObject> *bkgDataNode = new PHIODataNode<PHObject>(towerbackground, _backgroundName, "PHObject");
    bkgNode->addNode(bkgDataNode);
  }
  else
  {
    std::cout << PHWHERE << "::ERROR - " << _backgroundName << " pre-exists, but should not" << std::endl;
    exit(-1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void DetermineTowerBackground::FillNode(PHCompositeNode *topNode)
{

  TowerBackground *towerbackground = findNode::getClass<TowerBackground>(topNode, _backgroundName );
  if (!towerbackground)
  {
    std::cout << " ERROR -- can't find TowerBackground node after it should have been created" << std::endl;
    return;
  }
  else
  {

    towerbackground->set_UE( 0, _UE[0] );
    towerbackground->set_UE( 1, _UE[1] );
    towerbackground->set_UE( 2, _UE[2] );

    towerbackground->set_v2( _v2 );

    towerbackground->set_Psi2( _Psi2 );

  }


  return;
}
