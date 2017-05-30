#include "RetowerCEMC.h"

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

// standard includes
#include <iomanip>
#include <iostream>
#include <vector>

RetowerCEMC::RetowerCEMC(const std::string &name)
  : SubsysReco(name)
{


}

RetowerCEMC::~RetowerCEMC()
{
}

int RetowerCEMC::Init(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    std::cout << "RetowerCEMC::Init: initialized" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int RetowerCEMC::InitRun(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int RetowerCEMC::process_event(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    std::cout << "RetowerCEMC::process_event: entering" << std::endl;

  // pull out the tower containers and geometry objects at the start
  
  RawTowerContainer *towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  
  RawTowerGeomContainer *geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");

  if (verbosity > 0) {
    std::cout << "RetowerCEMC::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC towers" << std::endl;
  }

  // setup grid

  int NETA = geomIH->get_etabins();
  int NPHI = geomIH->get_phibins();
  
  _EMCAL_RETOWER_E.resize( NETA, std::vector<float>(NPHI, 0));

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

  // create new towers

  _emcal_retower = new RawTowerContainer( RawTowerDefs::CalorimeterId::HCALIN );

  for (int eta = 0; eta < NETA; eta++) {
    for (int phi = 0; phi < NPHI; phi++) {

      RawTower *new_tower = new RawTowerv1();

      new_tower->set_energy( _EMCAL_RETOWER_E[ eta ][ phi ] );
      _emcal_retower->AddTower( eta, phi, new_tower );
      
    }
  }

  // write out to node

  CreateNode(topNode);

  if (verbosity > 0) std::cout << "RetowerCEMC::process_event: exiting" << std::endl;

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
  
  PHIODataNode<PHObject> *emcalTowerNode = new PHIODataNode<PHObject>(_emcal_retower, "TOWER_CALIB_CEMC_RETOWER", "PHObject");
  emcalNode->addNode(emcalTowerNode);
  
  return Fun4AllReturnCodes::EVENT_OK;
}
