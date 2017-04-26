#include "CaloTriggerSim.h"

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
#include <g4cemc/RawTowerGeomContainer_Cylinderv1.h>

#include "CaloTriggerInfo_v1.h"

// standard includes
#include <iomanip>
#include <iostream>
#include <vector>

CaloTriggerSim::CaloTriggerSim(const std::string &name)
  : SubsysReco(name)
{
  // initiate sizes as -1 to tell module they can be set when it sees
  // the EMCal geometry for the first time

  _EMCAL_1x1_NETA = -1;
  _EMCAL_1x1_NPHI = -1;

  _EMCAL_2x2_NETA = -1;
  _EMCAL_2x2_NPHI = -1;

  _EMCAL_4x4_NETA = -1;
  _EMCAL_4x4_NPHI = -1;

  // these get cleared every event, but clear them at the
  // initiatlization step just in case

  _EMCAL_2x2_BEST_E = 0;
  _EMCAL_2x2_BEST_PHI = 0;
  _EMCAL_2x2_BEST_ETA = 0;

  _EMCAL_4x4_BEST_E = 0;
  _EMCAL_4x4_BEST_PHI = 0;
  _EMCAL_4x4_BEST_ETA = 0;
}

CaloTriggerSim::~CaloTriggerSim()
{
}

int CaloTriggerSim::Init(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    std::cout << "CaloTriggerSim::Init: initialized" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTriggerSim::InitRun(PHCompositeNode *topNode)
{
  return CreateNode(topNode);
}

int CaloTriggerSim::process_event(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    std::cout << "CaloTriggerSim::process_event: entering" << std::endl;

  RawTowerContainer *towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  if (verbosity > 0)
    std::cout << "CaloTriggerSim::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC towers" << std::endl;

  // get the binning from the geometry (different for 1D vs 2D...)
  RawTowerGeomContainer_Cylinderv1 *geomEM = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC");
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

    if (verbosity > 0)
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

    if (verbosity > 0 && tower->get_energy() > 1)
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

      // populate 2x2 map
      _EMCAL_2x2_MAP[ieta][iphi] = this_sum;

      // to calculate the eta, phi position, take the average of that of the 1x1's
      float this_eta = 0.5 * (geomEM->get_etacenter(2 * ieta) + geomEM->get_etacenter(2 * ieta + 1));
      float this_phi = 0.5 * (geomEM->get_phicenter(2 * iphi) + geomEM->get_phicenter(2 * iphi + 1));
      // wrap-around phi (apparently needed for 2D geometry?)
      if (this_phi > 3.14159) this_phi -= 2 * 3.14159;
      if (this_phi < -3.14159) this_phi += 2 * 3.14159;

      if (verbosity > 0 && this_sum > 1)
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

  if (verbosity > 0)
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
      this_sum += _EMCAL_2x2_MAP[ieta + 1][iphi];  // 2 * ieta + 1 is safe, since _EMCAL_4x4_NETA = _EMCAL_2x2_NETA - 1

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

      if (verbosity > 0 && this_sum > 1)
      {
        std::cout << "CaloTriggerSim::process_event: EMCal 4x4 tower eta ( bin ) / phi ( bin ) / E = " << std::setprecision(6) << this_eta << " ( " << ieta << " ) / " << this_phi << " ( " << iphi << " ) / " << this_sum << std::endl;
      }

      if (this_sum > _EMCAL_4x4_BEST_E)
      {
        _EMCAL_4x4_BEST_E = this_sum;
        _EMCAL_4x4_BEST_PHI = this_phi;
        _EMCAL_4x4_BEST_ETA = this_eta;
      }
    }
  }

  if (verbosity > 0)
  {
    std::cout << "CaloTriggerSim::process_event: best EMCal 4x4 window is at eta / phi = " << _EMCAL_4x4_BEST_ETA << " / " << _EMCAL_4x4_BEST_PHI << " and E = " << _EMCAL_4x4_BEST_E << std::endl;
  }

  FillNode(topNode);

  if (verbosity > 0) std::cout << "CaloTriggerSim::process_event: exiting" << std::endl;

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
  CaloTriggerInfo *triggerinfo = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");
  if (!triggerinfo)
  {
    triggerinfo = new CaloTriggerInfo_v1();
    PHIODataNode<PHObject> *TriggerNode = new PHIODataNode<PHObject>(triggerinfo, "CaloTriggerInfo", "PHObject");
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
  CaloTriggerInfo *triggerinfo = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");
  if (!triggerinfo)
  {
    std::cout << " ERROR -- can't find CaloTriggerInfo node after it should have been created" << std::endl;
    return;
  }
  else
  {
    triggerinfo->set_best_2x2_E(_EMCAL_2x2_BEST_E);
    triggerinfo->set_best_2x2_eta(_EMCAL_2x2_BEST_ETA);
    triggerinfo->set_best_2x2_phi(_EMCAL_2x2_BEST_PHI);

    triggerinfo->set_best_4x4_E(_EMCAL_4x4_BEST_E);
    triggerinfo->set_best_4x4_eta(_EMCAL_4x4_BEST_ETA);
    triggerinfo->set_best_4x4_phi(_EMCAL_4x4_BEST_PHI);
  }

  return;
}
