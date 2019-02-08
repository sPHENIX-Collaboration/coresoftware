#include "RetowerCEMC.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>

// sPHENIX includes
#include <calobase/RawTower.h>
#include <calobase/RawTowerv1.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>

// standard includes
#include <iomanip>
#include <iostream>
#include <vector>

RetowerCEMC::RetowerCEMC(const std::string &name)
  : SubsysReco(name), _NETA(-1), _NPHI(-1)
{


}

RetowerCEMC::~RetowerCEMC()
{
}

int RetowerCEMC::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    std::cout << "RetowerCEMC::Init: initialized" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int RetowerCEMC::InitRun(PHCompositeNode *topNode)
{

  CreateNode(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int RetowerCEMC::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    std::cout << "RetowerCEMC::process_event: entering" << std::endl;

  // pull out the tower containers and geometry objects at the start
  
  RawTowerContainer *towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");

  if (Verbosity() > 0) {
    std::cout << "RetowerCEMC::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC towers" << std::endl;
  }

  // setup grid

  if (_NETA < 0) {
    _NETA = geomIH->get_etabins();
    _NPHI = geomIH->get_phibins();

    _EMCAL_RETOWER_E.resize( _NETA, std::vector<float>(_NPHI, 0));
  }

  for (int eta = 0; eta < _NETA; eta++) {
    for (int phi = 0; phi < _NPHI; phi++) {
      _EMCAL_RETOWER_E[ eta ][ phi ] = 0;
    }
  }

  // partition existing CEMC energies among grid

  RawTowerContainer::ConstRange begin_end_EM = towersEM3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter) {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());
    
    int this_IHetabin = geomIH->get_etabin( tower_geom->get_eta() );
    int this_IHphibin = geomIH->get_phibin( tower_geom->get_phi() );
    float this_E = tower->get_energy();

    _EMCAL_RETOWER_E[ this_IHetabin ][ this_IHphibin ] += this_E;
    
  }

  RawTowerContainer* emcal_retower = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_CEMC_RETOWER");
  
  if (Verbosity() > 0) std::cout << "RetowerCEMC::process_event: filling TOWER_CALIB_CEMC_RETOWER node, with initial size = " << emcal_retower->size() << std::endl;

  // create new towers
  for (int eta = 0; eta < _NETA; eta++) {
    for (int phi = 0; phi < _NPHI; phi++) {

      RawTower *new_tower = new RawTowerv1();

      new_tower->set_energy( _EMCAL_RETOWER_E[ eta ][ phi ] );
      emcal_retower->AddTower( eta, phi, new_tower );
      
    }
  }

  if (Verbosity() > 0) std::cout << "RetowerCEMC::process_event: finished filling TOWER_CALIB_CEMC_RETOWER node, with final size = " << emcal_retower->size() << std::endl;

  if (Verbosity() > 0) std::cout << "RetowerCEMC::process_event: exiting" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int RetowerCEMC::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int RetowerCEMC::CreateNode(PHCompositeNode *topNode)
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

  RawTowerContainer* test_emcal_retower = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_CEMC_RETOWER");
  if ( !test_emcal_retower ) {

    if (Verbosity() > 0) std::cout << "RetowerCEMC::CreateNode : creating TOWER_CALIB_CEMC_RETOWER node " << std::endl;

    RawTowerContainer *emcal_retower = new RawTowerContainer( RawTowerDefs::CalorimeterId::HCALIN );
    PHIODataNode<PHObject> *emcalTowerNode = new PHIODataNode<PHObject>( emcal_retower, "TOWER_CALIB_CEMC_RETOWER", "PHObject");
    emcalNode->addNode(emcalTowerNode);
    
  } else {
    std::cout << "RetowerCEMC::CreateNode : TOWER_CALIB_CEMC_RETOWER already exists! " << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
