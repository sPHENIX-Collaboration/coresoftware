#include "SubtractTowersCS.h"

#include "TowerBackground.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/ConstituentSubtractor.hh>

// standard includes
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

SubtractTowersCS::SubtractTowersCS(const std::string &name)
  : SubsysReco(name)
  , _use_flow_modulation(false)
  , _alpha(1)
  , _DeltaRmax(0.3)
{
}

int SubtractTowersCS::InitRun(PHCompositeNode *topNode)
{
  CreateNode(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowersCS::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    std::cout << "SubtractTowersCS::process_event: entering, with _use_flow_modulation = " << _use_flow_modulation << std::endl;

  // pull out the tower containers and geometry objects at the start
  RawTowerContainer *towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
  RawTowerContainer *towersIH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");

  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowersCS::process_event: " << towersEM3->size() << " TOWER_CALIB_CEMC_RETOWER towers" << std::endl;
    std::cout << "SubtractTowersCS::process_event: " << towersIH3->size() << " TOWER_CALIB_HCALIN towers" << std::endl;
    std::cout << "SubtractTowersCS::process_event: " << towersOH3->size() << " TOWER_CALIB_HCALOUT towers" << std::endl;
  }

  // these should have already been created during InitRun()
  // note that these are the CS-versions
  RawTowerContainer *emcal_towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER_SUB1CS");
  RawTowerContainer *ihcal_towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN_SUB1CS");
  RawTowerContainer *ohcal_towers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT_SUB1CS");

  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowersCS::process_event: starting with " << emcal_towers->size() << " TOWER_CALIB_CEMC_RETOWER_SUB1CS towers" << std::endl;
    std::cout << "SubtractTowersCS::process_event: starting with " << ihcal_towers->size() << " TOWER_CALIB_HCALIN_SUB1CS towers" << std::endl;
    std::cout << "SubtractTowersCS::process_event: starting with " << ohcal_towers->size() << " TOWER_CALIB_HCALOUT_SUB1CS towers" << std::endl;
  }

  TowerBackground *towerbackground = findNode::getClass<TowerBackground>(topNode, "TowerBackground_Sub2");

  // read these in to use, even if we don't use flow modulation in the subtraction
  float background_v2 = towerbackground->get_v2();
  float background_Psi2 = towerbackground->get_Psi2();

  // set up constituent subtraction
  fastjet::contrib::ConstituentSubtractor subtractor;

  // free parameter for the type of distance between particle i and
  // ghost k. There are two options: "deltaR" or "angle" which are
  // defined as deltaR=sqrt((y_i-y_k)^2+(phi_i-phi_k)^2) or Euclidean
  // angle between the momenta
  subtractor.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);

  // free parameter for the maximal allowed distance between particle i and ghost k
  subtractor.set_max_distance(_DeltaRmax);

  // free parameter for the distance measure (the exponent of particle
  // pt). The larger the parameter alpha, the more are favoured the
  // lower pt particles in the subtraction process
  subtractor.set_alpha(_alpha);

  // free parameter for the density of ghosts. The smaller, the better
  // - but also the computation is slower.
  subtractor.set_ghost_area(0.01);

  // EM layer first
  std::vector<fastjet::PseudoJet> backgroundProxies_EM;
  std::vector<fastjet::PseudoJet> fullEvent_EM;
  std::vector<fastjet::PseudoJet> *backgroundProxies_EM_remaining = new std::vector<fastjet::PseudoJet>;

  // create PseudoJet version of EM towers in unsubtracted event
  RawTowerContainer::ConstRange begin_end_EM = towersEM3->getTowers();

  for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;

    double this_eta = geomIH->get_tower_geometry(tower->get_key())->get_eta();
    double this_phi = geomIH->get_tower_geometry(tower->get_key())->get_phi();
    double this_E = tower->get_energy();

    double this_pz = this_E * tanh(this_eta);
    double this_pt = sqrt(pow(this_E, 2) - pow(this_pz, 2));
    double this_px = this_pt * cos(this_phi);
    double this_py = this_pt * sin(this_phi);

    fastjet::PseudoJet this_pj(this_px, this_py, this_pz, this_E);

    fullEvent_EM.push_back(this_pj);
  }

  // create all new towers
  for (int eta = 0; eta < geomIH->get_etabins(); eta++)
  {
    for (int phi = 0; phi < geomIH->get_phibins(); phi++)
    {
      RawTower *new_tower = new RawTowerv1();
      new_tower->set_energy(0);
      emcal_towers->AddTower(eta, phi, new_tower);
    }
  }

  // create PseudoJet version of estimated background in EM layer
  for (RawTowerContainer::ConstIterator rtiter = emcal_towers->getTowers().first; rtiter != emcal_towers->getTowers().second; ++rtiter)
  {
    RawTower *tower = rtiter->second;

    double this_eta = geomIH->get_tower_geometry(tower->get_key())->get_eta();
    double this_phi = geomIH->get_tower_geometry(tower->get_key())->get_phi();

    double UE = towerbackground->get_UE(0).at(tower->get_bineta());
    if (_use_flow_modulation)
    {
      UE = UE * (1 + 2 * background_v2 * cos(2 * (this_phi - background_Psi2)));
    }

    double this_pz = UE * tanh(this_eta);
    double this_pt = sqrt(pow(UE, 2) - pow(this_pz, 2));
    double this_px = this_pt * cos(this_phi);
    double this_py = this_pt * sin(this_phi);

    fastjet::PseudoJet this_pj(this_px, this_py, this_pz, UE);

    backgroundProxies_EM.push_back(this_pj);

    if (Verbosity() > 5)
      std::cout << " SubtractTowersCS::process_event : background tower EM estimate for eta / phi = " << tower->get_bineta() << " / " << tower->get_binphi() << ", UE = " << UE << std::endl;
  }

  // constituent subtraction
  std::vector<fastjet::PseudoJet> correctedEvent_EM = subtractor.do_subtraction(fullEvent_EM, backgroundProxies_EM, backgroundProxies_EM_remaining);

  if (Verbosity() > 0)
  {
    std::cout << " SubtractTowersCS::process_event : vector lengths fullEvent_EM = " << fullEvent_EM.size() << " , backgroundProxies_EM = " << backgroundProxies_EM.size() << " , correctedEvent_EM = " << correctedEvent_EM.size() << ", backgroundProxies_EM_remaining = " << backgroundProxies_EM_remaining->size() << std::endl;

    double E_0 = 0, E_1 = 0, E_2 = 0, E_3 = 0;

    for (unsigned int n = 0; n < fullEvent_EM.size(); n++) E_0 += fullEvent_EM.at(n).E();
    for (unsigned int n = 0; n < backgroundProxies_EM.size(); n++) E_1 += backgroundProxies_EM.at(n).E();
    for (unsigned int n = 0; n < correctedEvent_EM.size(); n++) E_2 += correctedEvent_EM.at(n).E();
    for (unsigned int n = 0; n < backgroundProxies_EM_remaining->size(); n++) E_3 += backgroundProxies_EM_remaining->at(n).E();

    std::cout << "SubtractTowersCS::process_event EM : full event E - background E = " << E_0 << " - " << E_1 << " = " << E_0 - E_1 << ", subtracted E - remaining bkg E = " << E_2 << " - " << E_3 << " = " << E_2 - E_3 << std::endl;
  }

  // load subtracted towers into grid

  for (unsigned int n = 0; n < correctedEvent_EM.size(); n++)
  {
    float this_eta = correctedEvent_EM.at(n).eta();
    float this_phi = correctedEvent_EM.at(n).phi();
    float this_E = correctedEvent_EM.at(n).E();

    // look up tower by eta / phi...

    int this_etabin = geomIH->get_etabin(this_eta);
    int this_phibin = geomIH->get_phibin(this_phi);
    RawTower *tower = emcal_towers->getTower(this_etabin, this_phibin);
    tower->set_energy(this_E);

    if (Verbosity() > 5 || this_E < 0)
      std::cout << " SubtractTowersCS::process_event : creating subtracted EM tower for eta / phi = " << this_eta << " / " << this_phi << " ( " << tower->get_binphi() << " , " << this_phibin << " ) , sub. E  = " << this_E << (this_E < 0 ? " -- WARNING: negative E" : "") << std::endl;
  }

  // IHCal layer
  std::vector<fastjet::PseudoJet> backgroundProxies_IH;
  std::vector<fastjet::PseudoJet> fullEvent_IH;
  std::vector<fastjet::PseudoJet> *backgroundProxies_IH_remaining = new std::vector<fastjet::PseudoJet>;

  // create PseudoJet version of IH towers in unsubtracted event
  RawTowerContainer::ConstRange begin_end_IH = towersIH3->getTowers();

  for (RawTowerContainer::ConstIterator rtiter = begin_end_IH.first; rtiter != begin_end_IH.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;

    double this_eta = geomIH->get_tower_geometry(tower->get_key())->get_eta();
    double this_phi = geomIH->get_tower_geometry(tower->get_key())->get_phi();
    double this_E = tower->get_energy();

    double this_pz = this_E * tanh(this_eta);
    double this_pt = sqrt(pow(this_E, 2) - pow(this_pz, 2));
    double this_px = this_pt * cos(this_phi);
    double this_py = this_pt * sin(this_phi);

    fastjet::PseudoJet this_pj(this_px, this_py, this_pz, this_E);

    fullEvent_IH.push_back(this_pj);
  }

  // create all new towers
  for (int eta = 0; eta < geomIH->get_etabins(); eta++)
  {
    for (int phi = 0; phi < geomIH->get_phibins(); phi++)
    {
      RawTower *new_tower = new RawTowerv1();
      new_tower->set_energy(0);
      ihcal_towers->AddTower(eta, phi, new_tower);
    }
  }

  // create PseudoJet version of estimated background in IH layer
  for (RawTowerContainer::ConstIterator rtiter = ihcal_towers->getTowers().first; rtiter != ihcal_towers->getTowers().second; ++rtiter)
  {
    RawTower *tower = rtiter->second;

    double this_eta = geomIH->get_tower_geometry(tower->get_key())->get_eta();
    double this_phi = geomIH->get_tower_geometry(tower->get_key())->get_phi();

    double UE = towerbackground->get_UE(1).at(tower->get_bineta());
    if (_use_flow_modulation)
    {
      UE = UE * (1 + 2 * background_v2 * cos(2 * (this_phi - background_Psi2)));
    }

    double this_pz = UE * tanh(this_eta);
    double this_pt = sqrt(pow(UE, 2) - pow(this_pz, 2));
    double this_px = this_pt * cos(this_phi);
    double this_py = this_pt * sin(this_phi);

    fastjet::PseudoJet this_pj(this_px, this_py, this_pz, UE);

    backgroundProxies_IH.push_back(this_pj);

    if (Verbosity() > 5)
      std::cout << " SubtractTowersCS::process_event : background tower IH estimate for eta / phi = " << tower->get_bineta() << " / " << tower->get_binphi() << ", UE = " << UE << std::endl;
  }

  // constituent subtraction
  std::vector<fastjet::PseudoJet> correctedEvent_IH = subtractor.do_subtraction(fullEvent_IH, backgroundProxies_IH, backgroundProxies_IH_remaining);

  if (Verbosity() > 0)
  {
    std::cout << " SubtractTowersCS::process_event : vector lengths fullEvent_IH = " << fullEvent_IH.size() << " , backgroundProxies_IH = " << backgroundProxies_IH.size() << " , correctedEvent_IH = " << correctedEvent_IH.size() << ", backgroundProxies_IH_remaining = " << backgroundProxies_IH_remaining->size() << std::endl;

    double E_0 = 0, E_1 = 0, E_2 = 0, E_3 = 0;

    for (unsigned int n = 0; n < fullEvent_IH.size(); n++) E_0 += fullEvent_IH.at(n).E();
    for (unsigned int n = 0; n < backgroundProxies_IH.size(); n++) E_1 += backgroundProxies_IH.at(n).E();
    for (unsigned int n = 0; n < correctedEvent_IH.size(); n++) E_2 += correctedEvent_IH.at(n).E();
    for (unsigned int n = 0; n < backgroundProxies_IH_remaining->size(); n++) E_3 += backgroundProxies_IH_remaining->at(n).E();

    std::cout << "SubtractTowersCS::process_event IH : full event E - background E = " << E_0 << " - " << E_1 << " = " << E_0 - E_1 << ", subtracted E - remaining bkg E = " << E_2 << " - " << E_3 << " = " << E_2 - E_3 << std::endl;
  }

  // load subtracted towers into grid

  for (unsigned int n = 0; n < correctedEvent_IH.size(); n++)
  {
    float this_eta = correctedEvent_IH.at(n).eta();
    float this_phi = correctedEvent_IH.at(n).phi();
    float this_E = correctedEvent_IH.at(n).E();

    // look up tower by eta / phi...

    int this_etabin = geomIH->get_etabin(this_eta);
    int this_phibin = geomIH->get_phibin(this_phi);
    RawTower *tower = ihcal_towers->getTower(this_etabin, this_phibin);
    tower->set_energy(this_E);

    if (Verbosity() > 5 || this_E < 0)
      std::cout << " SubtractTowersCS::process_event : creating subtracted IH tower for eta / phi = " << this_eta << " / " << this_phi << " ( " << tower->get_binphi() << " , " << this_phibin << " ) , sub. E  = " << this_E << (this_E < 0 ? " -- WARNING: negative E" : "") << std::endl;
  }

  // OHCal layer
  std::vector<fastjet::PseudoJet> backgroundProxies_OH;
  std::vector<fastjet::PseudoJet> fullEvent_OH;
  std::vector<fastjet::PseudoJet> *backgroundProxies_OH_remaining = new std::vector<fastjet::PseudoJet>;

  // create PseudoJet version of OH towers in unsubtracted event
  RawTowerContainer::ConstRange begin_end_OH = towersOH3->getTowers();

  for (RawTowerContainer::ConstIterator rtiter = begin_end_OH.first; rtiter != begin_end_OH.second; ++rtiter)
  {
    RawTower *tower = rtiter->second;

    double this_eta = geomOH->get_tower_geometry(tower->get_key())->get_eta();
    double this_phi = geomOH->get_tower_geometry(tower->get_key())->get_phi();
    double this_E = tower->get_energy();

    double this_pz = this_E * tanh(this_eta);
    double this_pt = sqrt(pow(this_E, 2) - pow(this_pz, 2));
    double this_px = this_pt * cos(this_phi);
    double this_py = this_pt * sin(this_phi);

    fastjet::PseudoJet this_pj(this_px, this_py, this_pz, this_E);

    fullEvent_OH.push_back(this_pj);
  }

  // create all new towers
  for (int eta = 0; eta < geomOH->get_etabins(); eta++)
  {
    for (int phi = 0; phi < geomOH->get_phibins(); phi++)
    {
      RawTower *new_tower = new RawTowerv1();
      new_tower->set_energy(0);
      ohcal_towers->AddTower(eta, phi, new_tower);
    }
  }

  // create PseudoJet version of estimated background in OH layer
  for (RawTowerContainer::ConstIterator rtiter = ohcal_towers->getTowers().first; rtiter != ohcal_towers->getTowers().second; ++rtiter)
  {
    RawTower *tower = rtiter->second;

    double this_eta = geomOH->get_tower_geometry(tower->get_key())->get_eta();
    double this_phi = geomOH->get_tower_geometry(tower->get_key())->get_phi();

    double UE = towerbackground->get_UE(2).at(tower->get_bineta());
    if (_use_flow_modulation)
    {
      UE = UE * (1 + 2 * background_v2 * cos(2 * (this_phi - background_Psi2)));
    }

    double this_pz = UE * tanh(this_eta);
    double this_pt = sqrt(pow(UE, 2) - pow(this_pz, 2));
    double this_px = this_pt * cos(this_phi);
    double this_py = this_pt * sin(this_phi);

    fastjet::PseudoJet this_pj(this_px, this_py, this_pz, UE);

    backgroundProxies_OH.push_back(this_pj);

    if (Verbosity() > 5)
      std::cout << " SubtractTowersCS::process_event : background tower OH estimate for eta / phi = " << tower->get_bineta() << " / " << tower->get_binphi() << ", UE = " << UE << std::endl;
  }

  // constituent subtraction
  std::vector<fastjet::PseudoJet> correctedEvent_OH = subtractor.do_subtraction(fullEvent_OH, backgroundProxies_OH, backgroundProxies_OH_remaining);

  if (Verbosity() > 0)
  {
    std::cout << " SubtractTowersCS::process_event : vector lengths fullEvent_OH = " << fullEvent_OH.size() << " , backgroundProxies_OH = " << backgroundProxies_OH.size() << " , correctedEvent_OH = " << correctedEvent_OH.size() << ", backgroundProxies_OH_remaining = " << backgroundProxies_OH_remaining->size() << std::endl;

    double E_0 = 0, E_1 = 0, E_2 = 0, E_3 = 0;

    for (unsigned int n = 0; n < fullEvent_OH.size(); n++) E_0 += fullEvent_OH.at(n).E();
    for (unsigned int n = 0; n < backgroundProxies_OH.size(); n++) E_1 += backgroundProxies_OH.at(n).E();
    for (unsigned int n = 0; n < correctedEvent_OH.size(); n++) E_2 += correctedEvent_OH.at(n).E();
    for (unsigned int n = 0; n < backgroundProxies_OH_remaining->size(); n++) E_3 += backgroundProxies_OH_remaining->at(n).E();

    std::cout << "SubtractTowersCS::process_event OH : full event E - background E = " << E_0 << " - " << E_1 << " = " << E_0 - E_1 << ", subtracted E - remaining bkg E = " << E_2 << " - " << E_3 << " = " << E_2 - E_3 << std::endl;
  }

  // load subtracted towers into grid

  for (unsigned int n = 0; n < correctedEvent_OH.size(); n++)
  {
    float this_eta = correctedEvent_OH.at(n).eta();
    float this_phi = correctedEvent_OH.at(n).phi();
    float this_E = correctedEvent_OH.at(n).E();

    // look up tower by eta / phi...

    int this_etabin = geomOH->get_etabin(this_eta);
    int this_phibin = geomOH->get_phibin(this_phi);
    RawTower *tower = ohcal_towers->getTower(this_etabin, this_phibin);
    tower->set_energy(this_E);

    if (Verbosity() > 5 || this_E < 0)
      std::cout << " SubtractTowersCS::process_event : creating subtracted OH tower for eta / phi = " << this_eta << " / " << this_phi << " ( " << tower->get_binphi() << " , " << this_phibin << " ) , sub. E  = " << this_E << (this_E < 0 ? " -- WARNING: negative E" : "") << std::endl;
  }

  // cleanup

  delete backgroundProxies_EM_remaining;
  delete backgroundProxies_IH_remaining;
  delete backgroundProxies_OH_remaining;

  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowersCS::process_event: ending with " << emcal_towers->size() << " TOWER_CALIB_CEMC_RETOWER_SUB1CS towers" << std::endl;

    float EM_old_E = 0;
    float EM_new_E = 0;
    {
      RawTowerContainer::ConstRange begin_end_EM = towersEM3->getTowers();
      for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
      {
        RawTower *tower = rtiter->second;
        EM_old_E += tower->get_energy();
      }
    }
    {
      RawTowerContainer::ConstRange begin_end_EM = emcal_towers->getTowers();
      for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
      {
        RawTower *tower = rtiter->second;
        EM_new_E += tower->get_energy();
      }
    }
    std::cout << "SubtractTowersCS::process_event: old / new total E in EM layer = " << EM_old_E << " / " << EM_new_E << std::endl;

    std::cout << "SubtractTowersCS::process_event: ending with " << ihcal_towers->size() << " TOWER_CALIB_HCALIN_SUB1CS towers" << std::endl;

    float IH_old_E = 0;
    float IH_new_E = 0;
    {
      RawTowerContainer::ConstRange begin_end_EM = towersIH3->getTowers();
      for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
      {
        RawTower *tower = rtiter->second;
        IH_old_E += tower->get_energy();
      }
    }
    {
      RawTowerContainer::ConstRange begin_end_EM = ihcal_towers->getTowers();
      for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
      {
        RawTower *tower = rtiter->second;
        IH_new_E += tower->get_energy();
      }
    }
    std::cout << "SubtractTowersCS::process_event: old / new total E in IH layer = " << IH_old_E << " / " << IH_new_E << std::endl;

    std::cout << "SubtractTowersCS::process_event: ending with " << ohcal_towers->size() << " TOWER_CALIB_HCALOUT_SUB1CS towers" << std::endl;

    float OH_old_E = 0;
    float OH_new_E = 0;
    {
      RawTowerContainer::ConstRange begin_end_EM = towersOH3->getTowers();
      for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
      {
        RawTower *tower = rtiter->second;
        OH_old_E += tower->get_energy();
      }
    }
    {
      RawTowerContainer::ConstRange begin_end_EM = ohcal_towers->getTowers();
      for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
      {
        RawTower *tower = rtiter->second;
        OH_new_E += tower->get_energy();
      }
    }
    std::cout << "SubtractTowersCS::process_event: old / new total E in OH layer = " << OH_old_E << " / " << OH_new_E << std::endl;
  }

  if (Verbosity() > 0) std::cout << "SubtractTowersCS::process_event: exiting" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowersCS::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // store the new EMCal towers
  PHCompositeNode *emcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "CEMC"));
  if (!emcalNode)
  {
    std::cout << PHWHERE << "EMCal Node note found, doing nothing." << std::endl;
  }

  RawTowerContainer *test_emcal_tower = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER_SUB1CS");
  if (!test_emcal_tower)
  {
    if (Verbosity() > 0) std::cout << "SubtractTowersCS::CreateNode : creating TOWER_CALIB_CEMC_RETOWER_SUB1CS node " << std::endl;

    RawTowerContainer *emcal_towers = new RawTowerContainer(RawTowerDefs::CalorimeterId::HCALIN);
    PHIODataNode<PHObject> *emcalTowerNode = new PHIODataNode<PHObject>(emcal_towers, "TOWER_CALIB_CEMC_RETOWER_SUB1CS", "PHObject");
    emcalNode->addNode(emcalTowerNode);
  }
  else
  {
    std::cout << "SubtractTowersCS::CreateNode : TOWER_CALIB_CEMC_RETOWER_SUB1CS already exists! " << std::endl;
  }

  // store the new IHCal towers
  PHCompositeNode *ihcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "HCALIN"));
  if (!ihcalNode)
  {
    std::cout << PHWHERE << "IHCal Node note found, doing nothing." << std::endl;
  }

  RawTowerContainer *test_ihcal_tower = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN_SUB1CS");
  if (!test_ihcal_tower)
  {
    if (Verbosity() > 0) std::cout << "SubtractTowersCS::CreateNode : creating TOWER_CALIB_HCALIN_SUB1CS node " << std::endl;

    RawTowerContainer *ihcal_towers = new RawTowerContainer(RawTowerDefs::CalorimeterId::HCALIN);
    PHIODataNode<PHObject> *ihcalTowerNode = new PHIODataNode<PHObject>(ihcal_towers, "TOWER_CALIB_HCALIN_SUB1CS", "PHObject");
    ihcalNode->addNode(ihcalTowerNode);
  }
  else
  {
    std::cout << "SubtractTowersCS::CreateNode : TOWER_CALIB_HCALIN_SUB1CS already exists! " << std::endl;
  }

  // store the new OHCal towers
  PHCompositeNode *ohcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "HCALOUT"));
  if (!ohcalNode)
  {
    std::cout << PHWHERE << "OHCal Node note found, doing nothing." << std::endl;
  }

  RawTowerContainer *test_ohcal_tower = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT_SUB1CS");
  if (!test_ohcal_tower)
  {
    if (Verbosity() > 0) std::cout << "SubtractTowersCS::CreateNode : creating TOWER_CALIB_HCALOUT_SUB1CS node " << std::endl;

    RawTowerContainer *ohcal_towers = new RawTowerContainer(RawTowerDefs::CalorimeterId::HCALOUT);
    PHIODataNode<PHObject> *ohcalTowerNode = new PHIODataNode<PHObject>(ohcal_towers, "TOWER_CALIB_HCALOUT_SUB1CS", "PHObject");
    ohcalNode->addNode(ohcalTowerNode);
  }
  else
  {
    std::cout << "SubtractTowersCS::CreateNode : TOWER_CALIB_HCALOUT_SUB1CS already exists! " << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
