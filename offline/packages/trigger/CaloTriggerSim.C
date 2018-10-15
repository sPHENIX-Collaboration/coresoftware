#include "CaloTriggerSim.h"

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

#include "CaloTriggerInfo_v1.h"

// standard includes
#include <iomanip>
#include <iostream>
#include <vector>

CaloTriggerSim::CaloTriggerSim(const std::string &name)
  : SubsysReco(name), _emulate_truncation(0)
{
  // initiate sizes as -1 to tell module they can be set when it sees
  // the EMCal geometry for the first time

  _EMCAL_1x1_NETA = -1;
  _EMCAL_1x1_NPHI = -1;

  _EMCAL_2x2_NETA = -1;
  _EMCAL_2x2_NPHI = -1;

  _EMCAL_4x4_NETA = -1;
  _EMCAL_4x4_NPHI = -1;

  // do the same for full calo
  _FULLCALO_PHI_START = 0;
  _FULLCALO_PHI_END = 2*3.14159;
  
  _FULLCALO_0p1x0p1_NETA = -1;
  _FULLCALO_0p1x0p1_NPHI = -1;

  _FULLCALO_0p2x0p2_NETA = -1;
  _FULLCALO_0p2x0p2_NPHI = -1;

  _FULLCALO_0p4x0p4_NETA = -1;
  _FULLCALO_0p4x0p4_NPHI = -1;

  _FULLCALO_0p6x0p6_NETA = -1;
  _FULLCALO_0p6x0p6_NPHI = -1;

  _FULLCALO_0p8x0p8_NETA = -1;
  _FULLCALO_0p8x0p8_NPHI = -1;

  _FULLCALO_1p0x1p0_NETA = -1;
  _FULLCALO_1p0x1p0_NPHI = -1;

  // these get cleared every event, but clear them at the
  // initiatlization step just in case

  _EMCAL_2x2_BEST_E = 0;
  _EMCAL_2x2_BEST_PHI = 0;
  _EMCAL_2x2_BEST_ETA = 0;

  _EMCAL_4x4_BEST_E = 0;
  _EMCAL_4x4_BEST_PHI = 0;
  _EMCAL_4x4_BEST_ETA = 0;

  _EMCAL_4x4_BEST2_E = 0;
  _EMCAL_4x4_BEST2_PHI = 0;
  _EMCAL_4x4_BEST2_ETA = 0;

  _FULLCALO_0p2x0p2_BEST_E = 0;
  _FULLCALO_0p2x0p2_BEST_PHI = 0;
  _FULLCALO_0p2x0p2_BEST_ETA = 0;

  _FULLCALO_0p4x0p4_BEST_E = 0;
  _FULLCALO_0p4x0p4_BEST_PHI = 0;
  _FULLCALO_0p4x0p4_BEST_ETA = 0;

  _FULLCALO_0p6x0p6_BEST_E = 0;
  _FULLCALO_0p6x0p6_BEST_PHI = 0;
  _FULLCALO_0p6x0p6_BEST_ETA = 0;

  _FULLCALO_0p8x0p8_BEST_E = 0;
  _FULLCALO_0p8x0p8_BEST_PHI = 0;
  _FULLCALO_0p8x0p8_BEST_ETA = 0;

  _FULLCALO_1p0x1p0_BEST_E = 0;
  _FULLCALO_1p0x1p0_BEST_PHI = 0;
  _FULLCALO_1p0x1p0_BEST_ETA = 0;

}

CaloTriggerSim::~CaloTriggerSim()
{
}

void CaloTriggerSim::set_truncation( int emulate_truncation ) {
  _emulate_truncation = emulate_truncation;
}

float CaloTriggerSim::truncate_8bit( float raw_E  ) {

  if ( raw_E > 45.0 ) raw_E = 45.0;
  int counts = std::floor( raw_E / ( 45.0 / 256 ) );

  return counts * ( 45.0 / 256 );

}

int CaloTriggerSim::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    std::cout << "CaloTriggerSim::Init: initialized" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTriggerSim::InitRun(PHCompositeNode *topNode)
{
  return CreateNode(topNode);
}

int CaloTriggerSim::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    std::cout << "CaloTriggerSim::process_event: entering" << std::endl;

  // pull out the tower containers and geometry objects at the start
  RawTowerContainer *towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  RawTowerContainer *towersIH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");
  if (Verbosity() > 0) {
    std::cout << "CaloTriggerSim::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC towers" << std::endl;
    std::cout << "CaloTriggerSim::process_event: " << towersIH3->size() << " TOWER_CALIB_HCALIN towers" << std::endl;
    std::cout << "CaloTriggerSim::process_event: " << towersOH3->size() << " TOWER_CALIB_HCALOUT towers" << std::endl;
  }

  RawTowerGeomContainer_Cylinderv1 *geomEM = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC");
  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  // get the binning from the geometry (different for 1D vs 2D...)
  int geom_etabins = geomEM->get_etabins();
  int geom_phibins = geomEM->get_phibins();

  // if internal knowledge of geometry is unset, set it now (should
  // only happen once, on the first event)
  if (_EMCAL_1x1_NETA < 0)
  {
    _EMCAL_1x1_NETA = geom_etabins;
    _EMCAL_1x1_NPHI = geom_phibins;

    // half as many 2x2 windows along each axis as 1x1
    _EMCAL_2x2_NETA = geom_etabins / 2;
    _EMCAL_2x2_NPHI = geom_phibins / 2;

    // each 2x2 window defines a 4x4 window for which that 2x2 window
    // is the upper-left corner, so there are as many 4x4's as 2x2's
    // (except in eta, where the edge effect means there is 1 fewer)
    _EMCAL_4x4_NETA = geom_etabins / 2 - 1;
    _EMCAL_4x4_NPHI = geom_phibins / 2;

    // reset all maps
    _EMCAL_1x1_MAP.resize(_EMCAL_1x1_NETA, std::vector<float>(_EMCAL_1x1_NPHI, 0));
    _EMCAL_2x2_MAP.resize(_EMCAL_2x2_NETA, std::vector<float>(_EMCAL_2x2_NPHI, 0));
    _EMCAL_4x4_MAP.resize(_EMCAL_4x4_NETA, std::vector<float>(_EMCAL_4x4_NPHI, 0));

    if (Verbosity() > 0)
    {
      std::cout << "CaloTriggerSim::process_event: setting number of window in eta / phi,";
      std::cout << "1x1 are " << _EMCAL_1x1_NETA << " / " << _EMCAL_1x1_NPHI << ", ";
      std::cout << "2x2 are " << _EMCAL_2x2_NETA << " / " << _EMCAL_2x2_NPHI << ", ";
      std::cout << "4x4 are " << _EMCAL_4x4_NETA << " / " << _EMCAL_4x4_NPHI << std::endl;
    }
  }

  // reset 1x1 map
  for (int ieta = 0; ieta < _EMCAL_1x1_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _EMCAL_1x1_NPHI; iphi++)
    {
      _EMCAL_1x1_MAP[ieta][iphi] = 0;
    }
  }

  // iterate over EMCal towers, constructing 1x1's
  RawTowerContainer::ConstRange begin_end = towersEM3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());

    float this_eta = tower_geom->get_eta();
    float this_phi = tower_geom->get_phi();
    int this_etabin = geomEM->get_etabin(this_eta);
    int this_phibin = geomEM->get_phibin(this_phi);
    float this_E = tower->get_energy();

    _EMCAL_1x1_MAP[this_etabin][this_phibin] += this_E;

    if (Verbosity() > 1 && tower->get_energy() > 1)
    {
      std::cout << "CaloTriggerSim::process_event: EMCal 1x1 tower eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << this_etabin << " ) / " << this_phi << " ( " << this_phibin << " ) / " << this_E << std::endl;
    }
  }

  // reset 2x2 map and best
  for (int ieta = 0; ieta < _EMCAL_2x2_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _EMCAL_2x2_NPHI; iphi++)
    {
      _EMCAL_2x2_MAP[ieta][iphi] = 0;
    }
  }

  _EMCAL_2x2_BEST_E = 0;
  _EMCAL_2x2_BEST_PHI = 0;
  _EMCAL_2x2_BEST_ETA = 0;

  // now reconstruct 2x2 map from 1x1 map
  for (int ieta = 0; ieta < _EMCAL_2x2_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _EMCAL_2x2_NPHI; iphi++)
    {
      float this_sum = 0;

      this_sum += _EMCAL_1x1_MAP[2 * ieta][2 * iphi];
      this_sum += _EMCAL_1x1_MAP[2 * ieta][2 * iphi + 1];  // 2 * iphi + 1 is safe, since _EMCAL_2x2_NPHI = _EMCAL_1x1_NPHI / 2
      this_sum += _EMCAL_1x1_MAP[2 * ieta + 1][2 * iphi];  // 2 * ieta + 1 is safe, since _EMCAL_2x2_NETA = _EMCAL_1x1_NETA / 2
      this_sum += _EMCAL_1x1_MAP[2 * ieta + 1][2 * iphi + 1];

      if (_emulate_truncation) {
	this_sum = truncate_8bit( this_sum );
      }

      // populate 2x2 map
      _EMCAL_2x2_MAP[ieta][iphi] = this_sum;

      // to calculate the eta, phi position, take the average of that of the 1x1's
      float this_eta = 0.5 * (geomEM->get_etacenter(2 * ieta) + geomEM->get_etacenter(2 * ieta + 1));
      float this_phi = 0.5 * (geomEM->get_phicenter(2 * iphi) + geomEM->get_phicenter(2 * iphi + 1));
      // wrap-around phi (apparently needed for 2D geometry?)
      if (this_phi > 3.14159) this_phi -= 2 * 3.14159;
      if (this_phi < -3.14159) this_phi += 2 * 3.14159;

      if (Verbosity() > 1 && this_sum > 1)
      {
        std::cout << "CaloTriggerSim::process_event: EMCal 2x2 tower eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << ieta << " ) / " << this_phi << " ( " << iphi << " ) / " << this_sum << std::endl;
      }

      if (this_sum > _EMCAL_2x2_BEST_E)
      {
        _EMCAL_2x2_BEST_E = this_sum;
        _EMCAL_2x2_BEST_PHI = this_phi;
        _EMCAL_2x2_BEST_ETA = this_eta;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "CaloTriggerSim::process_event: best EMCal 2x2 window is at eta / phi = " << _EMCAL_2x2_BEST_ETA << " / " << _EMCAL_2x2_BEST_PHI << " and E = " << _EMCAL_2x2_BEST_E << std::endl;
  }

  // reset 4x4 map & best
  for (int ieta = 0; ieta < _EMCAL_4x4_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _EMCAL_4x4_NPHI; iphi++)
    {
      _EMCAL_4x4_MAP[ieta][iphi] = 0;
    }
  }

  _EMCAL_4x4_BEST_E = 0;
  _EMCAL_4x4_BEST_PHI = 0;
  _EMCAL_4x4_BEST_ETA = 0;

  int emcal_4x4_best_iphi = -1;
  int emcal_4x4_best_ieta = -1;

  // now reconstruct (sliding) 4x4 map from 2x2 map
  for (int ieta = 0; ieta < _EMCAL_4x4_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _EMCAL_4x4_NPHI; iphi++)
    {
      // for eta calculation (since eta distribution is potentially
      // non-uniform), average positions of all four towers
      float this_eta = 0.25 * (geomEM->get_etacenter(2 * ieta) + geomEM->get_etacenter(2 * ieta + 1) + geomEM->get_etacenter(2 * ieta + 2) + geomEM->get_etacenter(2 * ieta + 3));
      // for phi calculation (since phi distribution is uniform), take
      // first tower and add 1.5 tower widths
      float this_phi = geomEM->get_phicenter(2 * iphi) + 1.5 * (geomEM->get_phicenter(2 * iphi + 1) - geomEM->get_phicenter(2 * iphi));
      // wrap-around phi (apparently needed for 2D geometry?)
      if (this_phi > 3.14159) this_phi -= 2 * 3.14159;
      if (this_phi < -3.14159) this_phi += 2 * 3.14159;

      float this_sum = 0;

      this_sum += _EMCAL_2x2_MAP[ieta][iphi];
      this_sum += _EMCAL_2x2_MAP[ieta + 1][iphi];  // ieta + 1 is safe, since _EMCAL_4x4_NETA = _EMCAL_2x2_NETA - 1

      if (iphi != _EMCAL_4x4_NPHI - 1)
      {
        // if we are not in the last phi row, can safely access 'iphi+1'
        this_sum += _EMCAL_2x2_MAP[ieta][iphi + 1];
        this_sum += _EMCAL_2x2_MAP[ieta + 1][iphi + 1];
      }
      else
      {
        // if we are in the last phi row, wrap back around to zero
        this_sum += _EMCAL_2x2_MAP[ieta][0];
        this_sum += _EMCAL_2x2_MAP[ieta + 1][0];
      }

      _EMCAL_4x4_MAP[ieta][iphi] = this_sum;

      if (Verbosity() > 1 && this_sum > 1)
      {
        std::cout << "CaloTriggerSim::process_event: EMCal 4x4 tower eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << ieta << " ) / " << this_phi << " ( " << iphi << " ) / " << this_sum << std::endl;
      }

      if (this_sum > _EMCAL_4x4_BEST_E) {
	
	_EMCAL_4x4_BEST_E = this_sum;
	_EMCAL_4x4_BEST_PHI = this_phi;
	_EMCAL_4x4_BEST_ETA = this_eta;
	
	emcal_4x4_best_iphi = iphi;
	emcal_4x4_best_ieta = ieta;
  
      }
    }
  }


  _EMCAL_4x4_BEST2_E = 0;
  _EMCAL_4x4_BEST2_PHI = 0;
  _EMCAL_4x4_BEST2_ETA = 0;

  // find second-largest 4x4 which is > 1 tower away...
  for (int ieta = 0; ieta < _EMCAL_4x4_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _EMCAL_4x4_NPHI; iphi++)
    {

      int deta = ieta - emcal_4x4_best_ieta;
      int dphi = ( iphi - emcal_4x4_best_iphi ) % _EMCAL_4x4_NPHI ;

      if ( abs( deta ) < 1.5 && abs( dphi ) < 1.5 ) 
	continue;

      float this_eta = 0.25 * (geomEM->get_etacenter(2 * ieta) + geomEM->get_etacenter(2 * ieta + 1) + geomEM->get_etacenter(2 * ieta + 2) + geomEM->get_etacenter(2 * ieta + 3));
      float this_phi = geomEM->get_phicenter(2 * iphi) + 1.5 * (geomEM->get_phicenter(2 * iphi + 1) - geomEM->get_phicenter(2 * iphi));

      if (this_phi > 3.14159) this_phi -= 2 * 3.14159;
      if (this_phi < -3.14159) this_phi += 2 * 3.14159;

      float this_sum = _EMCAL_4x4_MAP[ieta][iphi];
      
      if (this_sum > _EMCAL_4x4_BEST2_E) {
	
	_EMCAL_4x4_BEST2_E = this_sum;
	_EMCAL_4x4_BEST2_PHI = this_phi;
	_EMCAL_4x4_BEST2_ETA = this_eta;

      }

    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "CaloTriggerSim::process_event: best EMCal 4x4 window is at eta / phi = " << _EMCAL_4x4_BEST_ETA << " / " << _EMCAL_4x4_BEST_PHI << " and E = " << _EMCAL_4x4_BEST_E << std::endl;
    std::cout << "CaloTriggerSim::process_event: 2nd best EMCal 4x4 window is at eta / phi = " << _EMCAL_4x4_BEST2_ETA << " / " << _EMCAL_4x4_BEST2_PHI << " and E = " << _EMCAL_4x4_BEST2_E << std::endl;
  }

  // begin full calo sim

  // get the 0.1x0.1 binning from the OHCal geometry
  int geomOH_etabins = geomOH->get_etabins();
  int geomOH_phibins = geomOH->get_phibins();

  // if internal knowledge of geometry is unset, set it now
  if (_FULLCALO_0p1x0p1_NETA < 0)
  {
    _FULLCALO_PHI_START = geomOH->get_phibounds( 0 ).first;
    _FULLCALO_PHI_END = geomOH->get_phibounds( geomOH_phibins - 1 ).second;

    _FULLCALO_0p1x0p1_NETA = geomOH_etabins;
    _FULLCALO_0p1x0p1_NPHI = geomOH_phibins;

    // half as many 0.2x0.2 windows along each axis as 0.1x0.1
    _FULLCALO_0p2x0p2_NETA = geomOH_etabins / 2;
    _FULLCALO_0p2x0p2_NPHI = geomOH_phibins / 2;

    // each 0.2x0.2 window defines a 0.4x0.4 window for which that
    // 0.2x0.2 window is the upper-left corner, so there are as many
    // 0.4x0.4's as 0.2x0.2's (except in eta, where the edge effect
    // means there is 1 fewer)
    _FULLCALO_0p4x0p4_NETA = geomOH_etabins / 2 - 1;
    _FULLCALO_0p4x0p4_NPHI = geomOH_phibins / 2;

    // for 0.6x0.6 windows, the above logic applies, except that the
    // edge effect causes there to be 2 fewer less in eta
    _FULLCALO_0p6x0p6_NETA = geomOH_etabins / 2 - 2;
    _FULLCALO_0p6x0p6_NPHI = geomOH_phibins / 2;

    // for 0.8x0.8 windows, the above logic applies, except that the
    // edge effect causes there to be 3 fewer less in eta
    _FULLCALO_0p8x0p8_NETA = geomOH_etabins / 2 - 3;
    _FULLCALO_0p8x0p8_NPHI = geomOH_phibins / 2;

    // for 1.0x1.0 windows, the above logic applies, except that the
    // edge effect causes there to be 4 fewer less in eta
    _FULLCALO_1p0x1p0_NETA = geomOH_etabins / 2 - 4;
    _FULLCALO_1p0x1p0_NPHI = geomOH_phibins / 2;

    // reset all maps
    _FULLCALO_0p1x0p1_MAP.resize(_FULLCALO_0p1x0p1_NETA, std::vector<float>(_FULLCALO_0p1x0p1_NPHI, 0));
    _FULLCALO_0p2x0p2_MAP.resize(_FULLCALO_0p2x0p2_NETA, std::vector<float>(_FULLCALO_0p2x0p2_NPHI, 0));
    _FULLCALO_0p4x0p4_MAP.resize(_FULLCALO_0p4x0p4_NETA, std::vector<float>(_FULLCALO_0p4x0p4_NPHI, 0));
    _FULLCALO_0p6x0p6_MAP.resize(_FULLCALO_0p6x0p6_NETA, std::vector<float>(_FULLCALO_0p6x0p6_NPHI, 0));
    _FULLCALO_0p8x0p8_MAP.resize(_FULLCALO_0p8x0p8_NETA, std::vector<float>(_FULLCALO_0p8x0p8_NPHI, 0));
    _FULLCALO_1p0x1p0_MAP.resize(_FULLCALO_1p0x1p0_NETA, std::vector<float>(_FULLCALO_1p0x1p0_NPHI, 0));

    if (Verbosity() > 0)
    {
      std::cout << "CaloTriggerSim::process_event: determining phi range for 0.1x0.1 full calo map: " << _FULLCALO_PHI_START << " to " << _FULLCALO_PHI_END << std::endl;
      std::cout << "CaloTriggerSim::process_event: setting number of full calo window in eta / phi:" << std::endl;
      std::cout << "  0.1x0.1 are " << _FULLCALO_0p1x0p1_NETA << " / " << _FULLCALO_0p1x0p1_NPHI << ", ";
      std::cout << "0.2x0.2 are " << _FULLCALO_0p2x0p2_NETA << " / " << _FULLCALO_0p2x0p2_NPHI << ", ";
      std::cout << "0.4x0.4 are " << _FULLCALO_0p4x0p4_NETA << " / " << _FULLCALO_0p4x0p4_NPHI << ", ";
      std::cout << "0.6x0.6 are " << _FULLCALO_0p6x0p6_NETA << " / " << _FULLCALO_0p6x0p6_NPHI << ", ";
      std::cout << "0.8x0.8 are " << _FULLCALO_0p8x0p8_NETA << " / " << _FULLCALO_0p8x0p8_NPHI << ", ";
      std::cout << "1.0x1.0 are " << _FULLCALO_1p0x1p0_NETA << " / " << _FULLCALO_1p0x1p0_NPHI << std::endl;
    }
  }

  // reset 0.1x0.1 map
  for (int ieta = 0; ieta < _FULLCALO_0p1x0p1_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_0p1x0p1_NPHI; iphi++)
    {
      _FULLCALO_0p1x0p1_MAP[ieta][iphi] = 0;
    }
  }

  // iterate over EMCal towers, filling in the 0.1x0.1 region they contribute to
  RawTowerContainer::ConstRange begin_end_EM = towersEM3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = geomEM->get_tower_geometry(tower->get_key());

    float this_eta = tower_geom->get_eta();
    float this_phi = tower_geom->get_phi();
    if (this_phi < _FULLCALO_PHI_START) this_phi += 2*3.14159;
    if (this_phi > _FULLCALO_PHI_END) this_phi -= 2*3.14159;
 
    // note: look up eta/phi index based on OHCal geometry, since this
    // defines the 0.1x0.1 regions
    int this_etabin = geomOH->get_etabin( this_eta );
    int this_phibin = geomOH->get_phibin( this_phi );
    float this_E = tower->get_energy();

    _FULLCALO_0p1x0p1_MAP[ this_etabin ][ this_phibin ] += this_E;

    if (Verbosity() > 1 && tower->get_energy() > 1)
    {
      std::cout << "CaloTriggerSim::process_event: EMCal tower at eta / phi (added to fullcalo map with etabin / phibin ) / E = " << std::setprecision(6) << this_eta << " / " << this_phi << " ( " << this_etabin << " / " << this_phibin << " ) / " << this_E << std::endl;
    }
  }

  // iterate over IHCal towers, filling in the 0.1x0.1 region they contribute to
  RawTowerContainer::ConstRange begin_end_IH = towersIH3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_IH.first; rtiter != begin_end_IH.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = geomIH->get_tower_geometry(tower->get_key());

    float this_eta = tower_geom->get_eta();
    float this_phi = tower_geom->get_phi();
    if (this_phi < _FULLCALO_PHI_START) this_phi += 2*3.14159;
    if (this_phi > _FULLCALO_PHI_END) this_phi -= 2*3.14159;
 
    // note: look up eta/phi index based on OHCal geometry, even though I
    // think it is by construction the same as the IHCal geometry...
    int this_etabin = geomOH->get_etabin( this_eta );
    int this_phibin = geomOH->get_phibin( this_phi );
    float this_E = tower->get_energy();

    _FULLCALO_0p1x0p1_MAP[ this_etabin ][ this_phibin ] += this_E;

    if (Verbosity() > 1 && tower->get_energy() > 0.5)
    {
      std::cout << "CaloTriggerSim::process_event: IHCal tower at eta / phi (added to fullcalo map with etabin / phibin ) / E = " << std::setprecision(6) << this_eta << " / " << this_phi << " ( " << this_etabin << " / " << this_phibin << " ) / " << this_E << std::endl;
    }
  }

  // iterate over OHCal towers, filling in the 0.1x0.1 region they contribute to
  RawTowerContainer::ConstRange begin_end_OH = towersOH3->getTowers();
  for (RawTowerContainer::ConstIterator rtiter = begin_end_OH.first; rtiter != begin_end_OH.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;
    RawTowerGeom *tower_geom = geomOH->get_tower_geometry(tower->get_key());

    float this_eta = tower_geom->get_eta();
    float this_phi = tower_geom->get_phi();
    if (this_phi < _FULLCALO_PHI_START) this_phi += 2*3.14159;
    if (this_phi > _FULLCALO_PHI_END) this_phi -= 2*3.14159;
 
    // note: use the nominal eta/phi index, since the fullcalo 0.1x0.1
    // map is defined by the OHCal geometry itself
    int this_etabin = geomOH->get_etabin( this_eta );
    int this_phibin = geomOH->get_phibin( this_phi );
    float this_E = tower->get_energy();

    _FULLCALO_0p1x0p1_MAP[ this_etabin ][ this_phibin ] += this_E;

    if (Verbosity() > 1 && tower->get_energy() > 0.5)
    {
      std::cout << "CaloTriggerSim::process_event: OHCal tower at eta / phi (added to fullcalo map with etabin / phibin ) / E = " << std::setprecision(6) << this_eta << " / " << this_phi << " ( " << this_etabin << " / " << this_phibin << " ) / " << this_E << std::endl;
    }
  }

  // reset 0.2x0.2 map and best
  for (int ieta = 0; ieta < _FULLCALO_0p2x0p2_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_0p2x0p2_NPHI; iphi++)
    {
      _FULLCALO_0p2x0p2_MAP[ ieta ][ iphi ] = 0;
    }
  }

  _FULLCALO_0p2x0p2_BEST_E = 0;
  _FULLCALO_0p2x0p2_BEST_PHI = 0;
  _FULLCALO_0p2x0p2_BEST_ETA = 0;

  // now reconstruct (non-sliding) 0.2x0.2 map from 0.1x0.1 map
  for (int ieta = 0; ieta < _FULLCALO_0p2x0p2_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_0p2x0p2_NPHI; iphi++)
    {
      float this_sum = 0;

      this_sum += _FULLCALO_0p1x0p1_MAP[2 * ieta][2 * iphi];
      this_sum += _FULLCALO_0p1x0p1_MAP[2 * ieta][2 * iphi + 1];  // 2 * iphi + 1 is safe, since _FULLCALO_0p2x0p2_NPHI = _FULLCALO_0p1x0p1_NPHI / 2
      this_sum += _FULLCALO_0p1x0p1_MAP[2 * ieta + 1][2 * iphi];  // 2 * ieta + 1 is safe, since _FULLCALO_0p2x0p2_NETA = _FULLCALO_0p1x0p1_NETA / 2
      this_sum += _FULLCALO_0p1x0p1_MAP[2 * ieta + 1][2 * iphi + 1];

      // populate 0.2x0.2 map
      _FULLCALO_0p2x0p2_MAP[ieta][iphi] = this_sum;

      // to calculate the eta, phi position, take the average of that
      // of the contributing 0.1x0.1's (which are defined by the OHCal geometry)
      float this_eta = 0.5 * (geomOH->get_etacenter(2 * ieta) + geomOH->get_etacenter(2 * ieta + 1));
      float this_phi = 0.5 * (geomOH->get_phicenter(2 * iphi) + geomOH->get_phicenter(2 * iphi + 1));

      if (Verbosity() > 1 && this_sum > 1)
      {
        std::cout << "CaloTriggerSim::process_event: FullCalo 0.2x0.2 window eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << ieta << " ) / " << this_phi << " ( " << iphi << " ) / " << this_sum << std::endl;
      }

      if (this_sum > _FULLCALO_0p2x0p2_BEST_E)
      {
        _FULLCALO_0p2x0p2_BEST_E = this_sum;
        _FULLCALO_0p2x0p2_BEST_PHI = this_phi;
        _FULLCALO_0p2x0p2_BEST_ETA = this_eta;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "CaloTriggerSim::process_event: best FullCalo 0.2x0.2 window is at eta / phi = " << _FULLCALO_0p2x0p2_BEST_ETA << " / " << _FULLCALO_0p2x0p2_BEST_PHI << " and E = " << _FULLCALO_0p2x0p2_BEST_E << std::endl;
  }

  // reset fullcalo 0.4x0.4 map & best
  for (int ieta = 0; ieta < _FULLCALO_0p4x0p4_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_0p4x0p4_NPHI; iphi++)
    {
      _FULLCALO_0p4x0p4_MAP[ ieta ][ iphi ] = 0;
    }
  }

  _FULLCALO_0p4x0p4_BEST_E = 0;
  _FULLCALO_0p4x0p4_BEST_PHI = 0;
  _FULLCALO_0p4x0p4_BEST_ETA = 0;

  // now reconstruct (sliding) 0.4x0.4 map from 0.2x0.2 map
  for (int ieta = 0; ieta < _FULLCALO_0p4x0p4_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_0p4x0p4_NPHI; iphi++)
    {
      // for eta calculation, use position of corner tower and add 1.5
      // tower widths
      float this_eta = geomOH->get_etacenter(2 * ieta) + 1.5 * ( geomOH->get_etacenter( 1 ) - geomOH->get_etacenter( 0 ) );
      // for phi calculation, use position of corner tower and add 1.5
      // tower widths
      float this_phi = geomOH->get_phicenter(2 * iphi) + 1.5 * (geomOH->get_phicenter( 1 ) - geomOH->get_phicenter( 0 ) );

      float this_sum = 0;

      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][iphi];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][iphi];  // 2 * ieta + 1 is safe, since _FULLCALO_0p4x0p4_NETA = _FULLCALO_0p4x0p4_NETA - 1

      // add 1 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];

      _FULLCALO_0p4x0p4_MAP[ieta][iphi] = this_sum;

      if (Verbosity() > 1 && this_sum > 2)
      {
        std::cout << "CaloTriggerSim::process_event: FullCalo  0.4x0.4 tower eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << ieta << " ) / " << this_phi << " ( " << iphi << " ) / " << this_sum << std::endl;
      }

      if (this_sum > _FULLCALO_0p4x0p4_BEST_E)
      {
        _FULLCALO_0p4x0p4_BEST_E = this_sum;
        _FULLCALO_0p4x0p4_BEST_PHI = this_phi;
        _FULLCALO_0p4x0p4_BEST_ETA = this_eta;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "CaloTriggerSim::process_event: best FullCalo 0.4x0.4 window is at eta / phi = " << _FULLCALO_0p4x0p4_BEST_ETA << " / " << _FULLCALO_0p4x0p4_BEST_PHI << " and E = " << _FULLCALO_0p4x0p4_BEST_E << std::endl;
  }

  // reset fullcalo 0.6x0.6 map & best
  for (int ieta = 0; ieta < _FULLCALO_0p6x0p6_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_0p6x0p6_NPHI; iphi++)
    {
      _FULLCALO_0p6x0p6_MAP[ ieta ][ iphi ] = 0;
    }
  }

  _FULLCALO_0p6x0p6_BEST_E = 0;
  _FULLCALO_0p6x0p6_BEST_PHI = 0;
  _FULLCALO_0p6x0p6_BEST_ETA = 0;

  // now reconstruct (sliding) 0.6x0.6 map from 0.2x0.2 map
  for (int ieta = 0; ieta < _FULLCALO_0p6x0p6_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_0p6x0p6_NPHI; iphi++)
    {
      // for eta calculation, use position of corner tower and add 2.5
      // tower widths
      float this_eta = geomOH->get_etacenter(2 * ieta) + 2.5 * ( geomOH->get_etacenter( 1 ) - geomOH->get_etacenter( 0 ) );
      // for phi calculation, use position of corner tower and add 2.5
      // tower widths
      float this_phi = geomOH->get_phicenter(2 * iphi) + 2.5 * (geomOH->get_phicenter( 1 ) - geomOH->get_phicenter( 0 ) );

      float this_sum = 0;

      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][iphi];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][iphi];  // ieta + 1 is safe, since _FULLCALO_0p6x0p6_NETA = _FULLCALO_0p2x0p2_NETA - 2
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][iphi];  // ieta + 2 is safe, since _FULLCALO_0p6x0p6_NETA = _FULLCALO_0p2x0p2_NETA - 2

      // add 1 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      // add 2 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];

      _FULLCALO_0p6x0p6_MAP[ieta][iphi] = this_sum;

      if (Verbosity() > 1 && this_sum > 3)
      {
        std::cout << "CaloTriggerSim::process_event: FullCalo  0.6x0.6 tower eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << ieta << " ) / " << this_phi << " ( " << iphi << " ) / " << this_sum << std::endl;
      }

      if (this_sum > _FULLCALO_0p6x0p6_BEST_E)
      {
        _FULLCALO_0p6x0p6_BEST_E = this_sum;
        _FULLCALO_0p6x0p6_BEST_PHI = this_phi;
        _FULLCALO_0p6x0p6_BEST_ETA = this_eta;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "CaloTriggerSim::process_event: best FullCalo 0.6x0.6 window is at eta / phi = " << _FULLCALO_0p6x0p6_BEST_ETA << " / " << _FULLCALO_0p6x0p6_BEST_PHI << " and E = " << _FULLCALO_0p6x0p6_BEST_E << std::endl;
  }

  // reset fullcalo 0.8x0.8 map & best
  for (int ieta = 0; ieta < _FULLCALO_0p8x0p8_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_0p8x0p8_NPHI; iphi++)
    {
      _FULLCALO_0p8x0p8_MAP[ ieta ][ iphi ] = 0;
    }
  }

  _FULLCALO_0p8x0p8_BEST_E = 0;
  _FULLCALO_0p8x0p8_BEST_PHI = 0;
  _FULLCALO_0p8x0p8_BEST_ETA = 0;

  // now reconstruct (sliding) 0.8x0.8 map from 0.2x0.2 map
  for (int ieta = 0; ieta < _FULLCALO_0p8x0p8_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_0p8x0p8_NPHI; iphi++)
    {
      // for eta calculation, use position of corner tower and add 3.5
      // tower widths
      float this_eta = geomOH->get_etacenter(2 * ieta) + 3.5 * ( geomOH->get_etacenter( 1 ) - geomOH->get_etacenter( 0 ) );
      // for phi calculation, use position of corner tower and add 3.5
      // tower widths
      float this_phi = geomOH->get_phicenter(2 * iphi) + 3.5 * (geomOH->get_phicenter( 1 ) - geomOH->get_phicenter( 0 ) );

      float this_sum = 0;

      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][iphi];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][iphi];  // ieta + 1 is safe, since _FULLCALO_0p8x0p8_NETA = _FULLCALO_0p2x0p2_NETA - 3
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][iphi];  // ieta + 2 is safe, since _FULLCALO_0p8x0p8_NETA = _FULLCALO_0p2x0p2_NETA - 3
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 3][iphi];  // ieta + 3 is safe, since _FULLCALO_0p8x0p8_NETA = _FULLCALO_0p2x0p2_NETA - 3

      // add 1 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 3][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      // add 2 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 3][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      // add 3 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 3 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 3 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][ ( iphi + 3 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 3][ ( iphi + 3 ) % _FULLCALO_0p2x0p2_NPHI ];

      _FULLCALO_0p8x0p8_MAP[ieta][iphi] = this_sum;

      if (Verbosity() > 1 && this_sum > 4)
      {
        std::cout << "CaloTriggerSim::process_event: FullCalo  0.8x0.8 tower eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << ieta << " ) / " << this_phi << " ( " << iphi << " ) / " << this_sum << std::endl;
      }

      if (this_sum > _FULLCALO_0p8x0p8_BEST_E)
      {
        _FULLCALO_0p8x0p8_BEST_E = this_sum;
        _FULLCALO_0p8x0p8_BEST_PHI = this_phi;
        _FULLCALO_0p8x0p8_BEST_ETA = this_eta;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "CaloTriggerSim::process_event: best FullCalo 0.8x0.8 window is at eta / phi = " << _FULLCALO_0p8x0p8_BEST_ETA << " / " << _FULLCALO_0p8x0p8_BEST_PHI << " and E = " << _FULLCALO_0p8x0p8_BEST_E << std::endl;
  }

  // reset fullcalo 1.0x1.0 map & best
  for (int ieta = 0; ieta < _FULLCALO_1p0x1p0_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_1p0x1p0_NPHI; iphi++)
    {
      _FULLCALO_1p0x1p0_MAP[ ieta ][ iphi ] = 0;
    }
  }

  _FULLCALO_1p0x1p0_BEST_E = 0;
  _FULLCALO_1p0x1p0_BEST_PHI = 0;
  _FULLCALO_1p0x1p0_BEST_ETA = 0;

  // now reconstruct (sliding) 1.0x1.0 map from 0.2x0.2 map
  for (int ieta = 0; ieta < _FULLCALO_1p0x1p0_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _FULLCALO_1p0x1p0_NPHI; iphi++)
    {
      // for eta calculation, use position of corner tower and add 4.5
      // tower widths
      float this_eta = geomOH->get_etacenter(2 * ieta) + 4.5 * ( geomOH->get_etacenter( 1 ) - geomOH->get_etacenter( 0 ) );
      // for phi calculation, use position of corner tower and add 4.5
      // tower widths
      float this_phi = geomOH->get_phicenter(2 * iphi) + 4.5 * (geomOH->get_phicenter( 1 ) - geomOH->get_phicenter( 0 ) );

      float this_sum = 0;

      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][iphi];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][iphi];  // ieta + 1 is safe, since _FULLCALO_1p0x1p0_NETA = _FULLCALO_0p2x0p2_NETA - 4
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][iphi];  // ieta + 2 is safe, since _FULLCALO_1p0x1p0_NETA = _FULLCALO_0p2x0p2_NETA - 4
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 3][iphi];  // ieta + 3 is safe, since _FULLCALO_1p0x1p0_NETA = _FULLCALO_0p2x0p2_NETA - 4
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 4][iphi];  // ieta + 4 is safe, since _FULLCALO_1p0x1p0_NETA = _FULLCALO_0p2x0p2_NETA - 4

      // add 1 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 3][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 4][ ( iphi + 1 ) % _FULLCALO_0p2x0p2_NPHI ];
      // add 2 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 3][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 4][ ( iphi + 2 ) % _FULLCALO_0p2x0p2_NPHI ];
      // add 3 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 3 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 3 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][ ( iphi + 3 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 3][ ( iphi + 3 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 4][ ( iphi + 3 ) % _FULLCALO_0p2x0p2_NPHI ];
      // add 4 to phi, but take modulus w.r.t. _FULLCALO_0p2x0p2_NPHI
      // in case we have wrapped back around
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta][ ( iphi + 4 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 1][ ( iphi + 4 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 2][ ( iphi + 4 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 3][ ( iphi + 4 ) % _FULLCALO_0p2x0p2_NPHI ];
      this_sum += _FULLCALO_0p2x0p2_MAP[ieta + 4][ ( iphi + 4 ) % _FULLCALO_0p2x0p2_NPHI ];

      _FULLCALO_1p0x1p0_MAP[ieta][iphi] = this_sum;

      if (Verbosity() > 1 && this_sum > 5)
      {
        std::cout << "CaloTriggerSim::process_event: FullCalo  1.0x1.0 tower eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << ieta << " ) / " << this_phi << " ( " << iphi << " ) / " << this_sum << std::endl;
      }

      if (this_sum > _FULLCALO_1p0x1p0_BEST_E)
      {
        _FULLCALO_1p0x1p0_BEST_E = this_sum;
        _FULLCALO_1p0x1p0_BEST_PHI = this_phi;
        _FULLCALO_1p0x1p0_BEST_ETA = this_eta;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "CaloTriggerSim::process_event: best FullCalo 1.0x1.0 window is at eta / phi = " << _FULLCALO_1p0x1p0_BEST_ETA << " / " << _FULLCALO_1p0x1p0_BEST_PHI << " and E = " << _FULLCALO_1p0x1p0_BEST_E << std::endl;
  }



  FillNode(topNode);

  if (Verbosity() > 0) std::cout << "CaloTriggerSim::process_event: exiting" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTriggerSim::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTriggerSim::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the trigger stuff under a sub-node directory
  PHCompositeNode *trigNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "TRIG"));
  if (!trigNode)
  {
    trigNode = new PHCompositeNode("TRIG");
    dstNode->addNode(trigNode);
  }

  // create the CaloTriggerInfo
  CaloTriggerInfo *triggerinfo = findNode::getClass<CaloTriggerInfo>(topNode, !_emulate_truncation ? "CaloTriggerInfo" : "CaloTriggerInfo_Truncate" );
  if (!triggerinfo)
  {
    triggerinfo = new CaloTriggerInfo_v1();
    PHIODataNode<PHObject> *TriggerNode = new PHIODataNode<PHObject>(triggerinfo, !_emulate_truncation ? "CaloTriggerInfo" : "CaloTriggerInfo_Truncate", "PHObject");
    trigNode->addNode(TriggerNode);
  }
  else
  {
    std::cout << PHWHERE << "::ERROR - CaloTriggerInfo pre-exists, but should not" << std::endl;
    exit(-1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTriggerSim::FillNode(PHCompositeNode *topNode)
{
  CaloTriggerInfo *triggerinfo = findNode::getClass<CaloTriggerInfo>(topNode, !_emulate_truncation ? "CaloTriggerInfo" : "CaloTriggerInfo_Truncate" );
  if (!triggerinfo)
  {
    std::cout << " ERROR -- can't find CaloTriggerInfo node after it should have been created" << std::endl;
    return;
  }
  else
  {
    triggerinfo->set_best_EMCal_2x2_E(_EMCAL_2x2_BEST_E);
    triggerinfo->set_best_EMCal_2x2_eta(_EMCAL_2x2_BEST_ETA);
    triggerinfo->set_best_EMCal_2x2_phi(_EMCAL_2x2_BEST_PHI);

    triggerinfo->set_best_EMCal_4x4_E(_EMCAL_4x4_BEST_E);
    triggerinfo->set_best_EMCal_4x4_eta(_EMCAL_4x4_BEST_ETA);
    triggerinfo->set_best_EMCal_4x4_phi(_EMCAL_4x4_BEST_PHI);

    triggerinfo->set_best2_EMCal_4x4_E(_EMCAL_4x4_BEST2_E);
    triggerinfo->set_best2_EMCal_4x4_eta(_EMCAL_4x4_BEST2_ETA);
    triggerinfo->set_best2_EMCal_4x4_phi(_EMCAL_4x4_BEST2_PHI);

    triggerinfo->set_best_FullCalo_0p2x0p2_E(_FULLCALO_0p2x0p2_BEST_E);
    triggerinfo->set_best_FullCalo_0p2x0p2_eta(_FULLCALO_0p2x0p2_BEST_ETA);
    triggerinfo->set_best_FullCalo_0p2x0p2_phi(_FULLCALO_0p2x0p2_BEST_PHI);

    triggerinfo->set_best_FullCalo_0p4x0p4_E(_FULLCALO_0p4x0p4_BEST_E);
    triggerinfo->set_best_FullCalo_0p4x0p4_eta(_FULLCALO_0p4x0p4_BEST_ETA);
    triggerinfo->set_best_FullCalo_0p4x0p4_phi(_FULLCALO_0p4x0p4_BEST_PHI);

    triggerinfo->set_best_FullCalo_0p6x0p6_E(_FULLCALO_0p6x0p6_BEST_E);
    triggerinfo->set_best_FullCalo_0p6x0p6_eta(_FULLCALO_0p6x0p6_BEST_ETA);
    triggerinfo->set_best_FullCalo_0p6x0p6_phi(_FULLCALO_0p6x0p6_BEST_PHI);

    triggerinfo->set_best_FullCalo_0p8x0p8_E(_FULLCALO_0p8x0p8_BEST_E);
    triggerinfo->set_best_FullCalo_0p8x0p8_eta(_FULLCALO_0p8x0p8_BEST_ETA);
    triggerinfo->set_best_FullCalo_0p8x0p8_phi(_FULLCALO_0p8x0p8_BEST_PHI);

    triggerinfo->set_best_FullCalo_1p0x1p0_E(_FULLCALO_1p0x1p0_BEST_E);
    triggerinfo->set_best_FullCalo_1p0x1p0_eta(_FULLCALO_1p0x1p0_BEST_ETA);
    triggerinfo->set_best_FullCalo_1p0x1p0_phi(_FULLCALO_1p0x1p0_BEST_PHI);

  }

  return;
}
