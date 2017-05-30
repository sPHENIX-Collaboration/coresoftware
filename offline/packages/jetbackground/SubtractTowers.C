#include "SubtractTowers.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>

// sPHENIX includes
#include <g4cemc/RawTower.h>
#include <g4cemc/RawTowerv1.h>
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

SubtractTowers::SubtractTowers(const std::string &name)
  : SubsysReco(name), _emcal_towers(nullptr), _ihcal_towers(nullptr), _ohcal_towers(nullptr)
{
  _background_iteration = 0;
}

SubtractTowers::~SubtractTowers()
{
}

int SubtractTowers::Init(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    std::cout << "SubtractTowers::Init: initialized" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowers::InitRun(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowers::process_event(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    std::cout << "SubtractTowers::process_event: entering" << std::endl;

  // pull out the tower containers and geometry objects at the start
  RawTowerContainer *towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
  RawTowerContainer *towersIH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");

  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  if (verbosity > 0) {
    std::cout << "SubtractTowers::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC_RETOWER towers" << std::endl;
    std::cout << "SubtractTowers::process_event: " << towersIH3->size() << " TOWER_CALIB_HCALIN towers" << std::endl;
    std::cout << "SubtractTowers::process_event: " << towersOH3->size() << " TOWER_CALIB_HCALOUT towers" << std::endl;
  }
  
  _emcal_towers = new RawTowerContainer( RawTowerDefs::CalorimeterId::HCALIN );
  _ihcal_towers = new RawTowerContainer( RawTowerDefs::CalorimeterId::HCALIN );
  _ohcal_towers = new RawTowerContainer( RawTowerDefs::CalorimeterId::HCALOUT );

  TowerBackground* towerbackground = findNode::getClass<TowerBackground>(topNode,"TowerBackground_Sub1");

  // EMCal
  
  // replicate existing towers
  RawTowerContainer::ConstRange begin_end_EM = towersEM3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    
    int this_etabin = tower->get_bineta();
    int this_phibin = tower->get_binphi();
    float this_E = tower->get_energy();

    RawTower *new_tower = new RawTowerv1();
    new_tower->set_energy( this_E );
    _emcal_towers->AddTower( this_etabin, this_phibin, new_tower );
    
  }

  // now fill in additional towers with zero energy to fill out the full grid
  // but note: after retowering, all of these should already exist...
  for (int eta = 0; eta < geomIH->get_etabins(); eta++) {
    for (int phi = 0; phi < geomIH->get_phibins(); phi++) {
      if ( ! _emcal_towers->getTower( eta, phi ) ) {
	RawTower *new_tower = new RawTowerv1();
	new_tower->set_energy( 0 );
	_emcal_towers->AddTower( eta, phi, new_tower );
      }
    }
  }

  // update towers for background subtraction...
  for (RawTowerContainer::ConstIterator rtiter = _emcal_towers->getTowers().first; rtiter != _emcal_towers->getTowers().second; ++rtiter) {
    RawTower *tower = rtiter->second;
    float raw_energy = tower->get_energy();
    float UE = towerbackground->get_UE( 0 ).at( tower->get_bineta() );
    float new_energy = raw_energy - UE;
    tower->set_energy( new_energy );
    if (verbosity > 5) 
      std::cout << " SubtractTowers::process_event : EMCal tower at eta / phi = " << tower->get_bineta() << " / " << tower->get_binphi() << ", pre-sub / after-sub E = " << raw_energy << " / " << tower->get_energy() << std::endl;
  }

  // IHCal

  // replicate existing towers
  RawTowerContainer::ConstRange begin_end_IH = towersIH3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_IH.first; rtiter != begin_end_IH.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());
    
    int this_etabin = geomIH->get_etabin( tower_geom->get_eta() );
    int this_phibin = geomIH->get_phibin( tower_geom->get_phi() );
    float this_E = tower->get_energy();

    RawTower *new_tower = new RawTowerv1();
    new_tower->set_energy( this_E );
    _ihcal_towers->AddTower( this_etabin, this_phibin, new_tower );
    
  }

  // now fill in additional towers with zero energy to fill out the full grid
  for (int eta = 0; eta < geomIH->get_etabins(); eta++) {
    for (int phi = 0; phi < geomIH->get_phibins(); phi++) {
      if ( ! _ihcal_towers->getTower( eta, phi ) ) {
	RawTower *new_tower = new RawTowerv1();
	new_tower->set_energy( 0 );
	_ihcal_towers->AddTower( eta, phi, new_tower );
      }
    }
  }

  // update towers for background subtraction...
  for (RawTowerContainer::ConstIterator rtiter = _ihcal_towers->getTowers().first; rtiter != _ihcal_towers->getTowers().second; ++rtiter) {
    RawTower *tower = rtiter->second;
    float raw_energy = tower->get_energy();
    float UE = towerbackground->get_UE( 1 ).at( tower->get_bineta() );
    float new_energy = raw_energy - UE;
    tower->set_energy( new_energy );
  }

  // OHCal

  // replicate existing towers
  RawTowerContainer::ConstRange begin_end_OH = towersOH3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_OH.first; rtiter != begin_end_OH.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());
    
    int this_etabin = geomOH->get_etabin( tower_geom->get_eta() );
    int this_phibin = geomOH->get_phibin( tower_geom->get_phi() );
    float this_E = tower->get_energy();

    RawTower *new_tower = new RawTowerv1();
    new_tower->set_energy( this_E );
    _ohcal_towers->AddTower( this_etabin, this_phibin, new_tower );
    
  }

  // now fill in additional towers with zero energy to fill out the full grid
  for (int eta = 0; eta < geomOH->get_etabins(); eta++) {
    for (int phi = 0; phi < geomOH->get_phibins(); phi++) {
      if ( ! _ohcal_towers->getTower( eta, phi ) ) {
	RawTower *new_tower = new RawTowerv1();
	new_tower->set_energy( 0 );
	_ohcal_towers->AddTower( eta, phi, new_tower );
      }
    }
  }

  // update towers for background subtraction...
  for (RawTowerContainer::ConstIterator rtiter = _ohcal_towers->getTowers().first; rtiter != _ohcal_towers->getTowers().second; ++rtiter) {
    RawTower *tower = rtiter->second;
    float raw_energy = tower->get_energy();
    float UE = towerbackground->get_UE( 2 ).at( tower->get_bineta() );
    float new_energy = raw_energy - UE;
    tower->set_energy( new_energy );
  }

  if (verbosity > 0) {
    std::cout << "SubtractTowers::process_event: " << _emcal_towers->size() << " TOWER_CALIB_CEMC_RETOWER_SUB1 towers" << std::endl;
    std::cout << "SubtractTowers::process_event: " << _ihcal_towers->size() << " TOWER_CALIB_HCALIN_SUB1 towers" << std::endl;
    std::cout << "SubtractTowers::process_event: " << _ohcal_towers->size() << " TOWER_CALIB_HCALOUT_SUB1 towers" << std::endl;
  }

  CreateNode(topNode);

  if (verbosity > 0) std::cout << "SubtractTowers::process_event: exiting" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowers::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowers::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the new EMCal towers
  PHCompositeNode *emcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "CEMC"));
  if (!emcalNode) {
    std::cout << PHWHERE << "EMCal Node note found, doing nothing." << std::endl;
  }
  
  PHIODataNode<PHObject> *emcalTowerNode = new PHIODataNode<PHObject>(_emcal_towers, "TOWER_CALIB_CEMC_RETOWER_SUB1", "PHObject");
  emcalNode->addNode(emcalTowerNode);

  // store the new IHCal towers
  PHCompositeNode *ihcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "HCALIN"));
  if (!ihcalNode) {
    std::cout << PHWHERE << "IHCal Node note found, doing nothing." << std::endl;
  }
  
  PHIODataNode<PHObject> *ihcalTowerNode = new PHIODataNode<PHObject>(_ihcal_towers, "TOWER_CALIB_HCALIN_SUB1", "PHObject");
  ihcalNode->addNode(ihcalTowerNode);

  // store the new OHCal towers
  PHCompositeNode *ohcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "HCALOUT"));
  if (!ohcalNode) {
    std::cout << PHWHERE << "OHCal Node note found, doing nothing." << std::endl;
  }

  PHIODataNode<PHObject> *ohcalTowerNode = new PHIODataNode<PHObject>(_ohcal_towers, "TOWER_CALIB_HCALOUT_SUB1", "PHObject");
  ohcalNode->addNode(ohcalTowerNode);

  return Fun4AllReturnCodes::EVENT_OK;
}
