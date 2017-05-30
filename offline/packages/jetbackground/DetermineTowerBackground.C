#include "DetermineTowerBackground.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>

// sPHENIX includes
#include <g4cemc/RawTower.h>
#include <g4cemc/RawTowerContainer.h>
#include <g4cemc/RawTowerGeom.h>
#include <g4cemc/RawTowerGeomContainer.h>
#include <g4cemc/RawTowerGeomContainer_Cylinderv1.h>

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include "TowerBackground_v1.h"

// standard includes
#include <iomanip>
#include <iostream>
#include <vector>

DetermineTowerBackground::DetermineTowerBackground(const std::string &name)
  : SubsysReco(name)
{

  _v2[0] = 0;
  _v2[1] = 0;
  _v2[2] = 0;

  _Psi2[0] = 0;
  _Psi2[1] = 0;
  _Psi2[2] = 0;

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
  if (verbosity > 0)
    std::cout << "DetermineTowerBackground::Init: initialized" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int DetermineTowerBackground::InitRun(PHCompositeNode *topNode)
{
  return CreateNode(topNode);
}

int DetermineTowerBackground::process_event(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    std::cout << "DetermineTowerBackground::process_event: entering" << std::endl;

  // clear seed eta/phi positions
  _seed_eta.resize(0);
  _seed_phi.resize(0);

  // 
  if (_seed_type == 1) {
    JetMap* reco2_jets = findNode::getClass<JetMap>(topNode,"AntiKt_Tower_r02");

    if (verbosity > 1)
      std::cout << "DetermineTowerBackground::proess_event: examining possible seeds ... " << std::endl;
    
    for (JetMap::Iter iter = reco2_jets->begin(); iter != reco2_jets->end(); ++iter) {
      Jet* this_jet = iter->second;
      
      float this_pt = this_jet->get_pt();
      float this_phi = this_jet->get_phi();
      float this_eta = this_jet->get_eta();

      if (this_jet->get_pt() < 25) continue;

      _seed_eta.push_back( this_eta );
      _seed_phi.push_back( this_phi );

      if (verbosity > 1)
	std::cout << "DetermineTowerBackground::proess_event: adding seed at eta / phi = " << this_eta << " / " << this_phi << " ( R=0.2 jet with pt = " << this_pt << " ) " << std::endl;
    }
    
  }

  // pull out the tower containers and geometry objects at the start
  RawTowerContainer *towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
  RawTowerContainer *towersIH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");
  if (verbosity > 0) {
    std::cout << "DetermineTowerBackground::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC_RETOWER towers" << std::endl;
    std::cout << "DetermineTowerBackground::process_event: " << towersIH3->size() << " TOWER_CALIB_HCALIN towers" << std::endl;
    std::cout << "DetermineTowerBackground::process_event: " << towersOH3->size() << " TOWER_CALIB_HCALOUT towers" << std::endl;
  }

  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

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

    if (verbosity > 0) {
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
    
    if (verbosity > 1 && tower->get_energy() > 1)
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

    if (verbosity > 1 && tower->get_energy() > 1)
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

    if (verbosity > 1 && tower->get_energy() > 1)
    {
      std::cout << "DetermineTowerBackground::process_event: OHCal tower at eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << this_etabin << " ) / " << this_phi << " ( " << this_phibin << " ) / " << this_E << std::endl;
    }
  }

  // calculate energy densities...

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
	    if (verbosity > 10) std::cout << " setting excluded mark from seed at eta / phi = " << _seed_eta[ iseed ] << " / " << _seed_phi[ iseed ] << std::endl;
	  }
	}

	if (!isExcluded) {
	  if ( layer == 0 ) total_E += _EMCAL_E[ eta ][ phi ];
	  if ( layer == 1 ) total_E += _IHCAL_E[ eta ][ phi ];
	  if ( layer == 2 ) total_E += _OHCAL_E[ eta ][ phi ];
	  total_tower++;
	} else {
	  if (verbosity > 10) std::cout << " tower at eta / phi = " << this_eta << " / " << this_phi << " with E = " << total_E << " excluded due to seed " << std::endl;
	}
      }

      std::pair< float, float > etabounds = geomIH->get_etabounds( eta );
      std::pair< float, float > phibounds = geomIH->get_phibounds( 0 );

      float deta = etabounds.second - etabounds.first;
      float dphi = phibounds.second - phibounds.first;
      float total_area = total_tower * deta * dphi;
      _UE[ layer ].at( eta ) = total_E / total_tower;
      
      if (verbosity > 1 ) {
	std::cout << "DetermineTowerBackground::process_event: at layer / eta index ( eta range ) = " << layer << " / " << eta << " ( " << etabounds.first << " - " << etabounds.second << " ) , total E / total Ntower / total area = " << total_E << " / " << total_tower << " / " << total_area << " , UE per tower = " << total_E / total_tower << std::endl;
      }
      
    }
  }
  
  // 

  FillNode(topNode);

  if (verbosity > 0) std::cout << "DetermineTowerBackground::process_event: exiting" << std::endl;

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

    towerbackground->set_v2( 0, _v2[0] );
    towerbackground->set_v2( 0, _v2[1] );
    towerbackground->set_v2( 0, _v2[2] );

    towerbackground->set_Psi2( 0, _Psi2[0] );
    towerbackground->set_Psi2( 0, _Psi2[1] );
    towerbackground->set_Psi2( 0, _Psi2[2] );

  }


  return;
}
