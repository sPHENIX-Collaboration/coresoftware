#include "SubtractTowers.h"

#include "TowerBackground.h"

// sPHENIX includes
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// standard includes
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

SubtractTowers::SubtractTowers(const std::string &name)
  : SubsysReco(name)
{
}

int SubtractTowers::InitRun(PHCompositeNode *topNode)
{
  CreateNode(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowers::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowers::process_event: entering, with _use_flow_modulation = " << _use_flow_modulation << std::endl;
  }

  // pull out the tower containers and geometry objects at the start
  EMTowerName = emcal_input_node();
  IHTowerName = ihcal_input_node();
  OHTowerName = ohcal_input_node();
  TowerInfoContainer *towerinfosEM3 = findNode::getClass<TowerInfoContainer>(topNode, EMTowerName);
  TowerInfoContainer *towerinfosIH3 = findNode::getClass<TowerInfoContainer>(topNode, IHTowerName);
  TowerInfoContainer *towerinfosOH3 = findNode::getClass<TowerInfoContainer>(topNode, OHTowerName);
  if (!towerinfosEM3 || !towerinfosIH3 || !towerinfosOH3)
  {
    std::cout << PHWHERE << " missing input tower node " << EMTowerName << ", " << IHTowerName << " or " << OHTowerName << ", exiting" << std::endl;
    exit(1);
  }

  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowers::process_event: " << towerinfosEM3->size() << EMTowerName << " towers" << std::endl;
    std::cout << "SubtractTowers::process_event: " << towerinfosIH3->size() << IHTowerName << " towers" << std::endl;
    std::cout << "SubtractTowers::process_event: " << towerinfosOH3->size() << OHTowerName << " towers" << std::endl;
  }

  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, m_ihcal_geom_node);
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, m_ohcal_geom_node);
  if (_use_flow_modulation && (!geomIH || !geomOH))
  {
    std::cout << PHWHERE << " missing tower geometry " << m_ihcal_geom_node << " or " << m_ohcal_geom_node << ", exiting" << std::endl;
    exit(1);
  }

  // these should have already been created during InitRun()
  EMTowerName = emcal_output_node();
  IHTowerName = ihcal_output_node();
  OHTowerName = ohcal_output_node();
  TowerInfoContainer *emcal_towerinfos = findNode::getClass<TowerInfoContainer>(topNode, EMTowerName);
  TowerInfoContainer *ihcal_towerinfos = findNode::getClass<TowerInfoContainer>(topNode, IHTowerName);
  TowerInfoContainer *ohcal_towerinfos = findNode::getClass<TowerInfoContainer>(topNode, OHTowerName);
  if (!emcal_towerinfos || !ihcal_towerinfos || !ohcal_towerinfos)
  {
    std::cout << PHWHERE << " missing output tower node " << EMTowerName << ", " << IHTowerName << " or " << OHTowerName << ", exiting" << std::endl;
    exit(1);
  }

  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowers::process_event: starting with " << emcal_towerinfos->size() << EMTowerName << " towers" << std::endl;
    std::cout << "SubtractTowers::process_event: starting with " << ihcal_towerinfos->size() << IHTowerName << " towers" << std::endl;
    std::cout << "SubtractTowers::process_event: starting with " << ohcal_towerinfos->size() << OHTowerName << " towers" << std::endl;
  }

  TowerBackground *towerbackground = findNode::getClass<TowerBackground>(topNode, m_towerBackgroundNode);
  if (!towerbackground)
  {
    std::cout << PHWHERE << " TowerBackground node " << m_towerBackgroundNode << " not found (see set_inputTowerBackgroundNode), exiting" << std::endl;
    exit(1);
  }
  // read these in to use, even if we don't use flow modulation in the subtraction
  float background_v2 = towerbackground->get_v2();
  float background_Psi2 = towerbackground->get_Psi2();

  // EMCal

  // replicate existing towers
  unsigned int nchannels_em = towerinfosEM3->size();
  for (unsigned int channel = 0; channel < nchannels_em; channel++)
  {
    TowerInfo *tower = towerinfosEM3->get_tower_at_channel(channel);
    unsigned int towerkey = towerinfosEM3->encode_key(channel);
    int ieta = towerinfosEM3->getTowerEtaBin(towerkey);
    int iphi = towerinfosEM3->getTowerPhiBin(towerkey);
    float raw_energy = tower->get_energy();
    float UE = towerbackground->get_UE(0).at(ieta);
    if (_use_flow_modulation)
    {
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
      float tower_phi = geomIH->get_tower_geometry(key)->get_phi();
      float modulation_factor = 1 + 2 * background_v2 * std::cos(2 * (tower_phi - background_Psi2));
      modulation_factor = std::max(0.F, modulation_factor);
      UE = UE * modulation_factor;
    }
    float new_energy = raw_energy - UE;
    // if a tower is masked, leave it at zero
    if (!tower->get_isGood())
    {
      new_energy = 0;
    }

    // the SUB1 container is cloned once in CreateNode; without this its status
    // bits are whatever the clone source held at InitRun (event 0 only)
    emcal_towerinfos->get_tower_at_channel(channel)->set_status(tower->get_status());
    emcal_towerinfos->get_tower_at_channel(channel)->set_time(tower->get_time());
    emcal_towerinfos->get_tower_at_channel(channel)->set_energy(new_energy);

    if (Verbosity() > 5)
    {
      std::cout << " SubtractTowers::process_event : EMCal tower at ieta / iphi = " << ieta << " / " << iphi << ", pre-sub / after-sub E = " << raw_energy << " / " << new_energy << std::endl;
    }
  }

  // IHCal
  // replicate existing towers
  unsigned int nchannels_ih = towerinfosIH3->size();
  for (unsigned int channel = 0; channel < nchannels_ih; channel++)
  {
    TowerInfo *tower = towerinfosIH3->get_tower_at_channel(channel);
    unsigned int towerkey = towerinfosIH3->encode_key(channel);
    int ieta = towerinfosIH3->getTowerEtaBin(towerkey);
    int iphi = towerinfosIH3->getTowerPhiBin(towerkey);

    float raw_energy = tower->get_energy();
    float UE = towerbackground->get_UE(1).at(ieta);
    if (_use_flow_modulation)
    {
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, ieta, iphi);
      float tower_phi = geomIH->get_tower_geometry(key)->get_phi();
      float modulation_factor = 1 + 2 * background_v2 * std::cos(2 * (tower_phi - background_Psi2));
      modulation_factor = std::max(0.F, modulation_factor);
      UE = UE * modulation_factor;
    }
    float new_energy = raw_energy - UE;
    // if a tower is masked, leave it at zero
    if (!tower->get_isGood())
    {
      new_energy = 0;
    }

    // the SUB1 container is cloned once in CreateNode; without this its status
    // bits are whatever the clone source held at InitRun (event 0 only)
    ihcal_towerinfos->get_tower_at_channel(channel)->set_status(tower->get_status());
    ihcal_towerinfos->get_tower_at_channel(channel)->set_time(tower->get_time());
    ihcal_towerinfos->get_tower_at_channel(channel)->set_energy(new_energy);
    if (Verbosity() > 5)
    {
      std::cout << "SubtractTowers::process_event : IHCal tower at ieta / iphi = " << ieta << " / " << iphi << ", pre-sub / after-sub E = " << raw_energy << " / " << new_energy << std::endl;
    }
  }

  // OHCal
  // replicate existing towers
  unsigned int nchannels_oh = towerinfosOH3->size();
  for (unsigned int channel = 0; channel < nchannels_oh; channel++)
  {
    TowerInfo *tower = towerinfosOH3->get_tower_at_channel(channel);
    unsigned int towerkey = towerinfosOH3->encode_key(channel);
    int ieta = towerinfosOH3->getTowerEtaBin(towerkey);
    int iphi = towerinfosOH3->getTowerPhiBin(towerkey);
    float raw_energy = tower->get_energy();
    float UE = towerbackground->get_UE(2).at(ieta);
    if (_use_flow_modulation)
    {
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
      float tower_phi = geomOH->get_tower_geometry(key)->get_phi();
      float modulation_factor = 1 + 2 * background_v2 * std::cos(2 * (tower_phi - background_Psi2));
      modulation_factor = std::max(0.F, modulation_factor);
      UE = UE * modulation_factor;
    }
    float new_energy = raw_energy - UE;
    // if a tower is masked, leave it at zero
    if (!tower->get_isGood())
    {
      new_energy = 0;
    }

    // the SUB1 container is cloned once in CreateNode; without this its status
    // bits are whatever the clone source held at InitRun (event 0 only)
    ohcal_towerinfos->get_tower_at_channel(channel)->set_status(tower->get_status());
    ohcal_towerinfos->get_tower_at_channel(channel)->set_time(tower->get_time());
    ohcal_towerinfos->get_tower_at_channel(channel)->set_energy(new_energy);
    if (Verbosity() > 5)
    {
      std::cout << "SubtractTowers::process_event : OHCal tower at ieta / iphi = " << ieta << " / " << iphi << ", pre-sub / after-sub E = " << raw_energy << " / " << new_energy << std::endl;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowers::process_event: ending with " << emcal_towerinfos->size() << " " << emcal_output_node() << " towers" << std::endl;
    std::cout << "SubtractTowers::process_event: ending with " << ihcal_towerinfos->size() << " " << ihcal_output_node() << " towers" << std::endl;
    std::cout << "SubtractTowers::process_event: ending with " << ohcal_towerinfos->size() << " " << ohcal_output_node() << " towers" << std::endl;
  }

  if (Verbosity() > 0)
  {
    std::cout << "SubtractTowers::process_event: exiting" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SubtractTowers::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  IHTowerName = ihcal_input_node();
  TowerInfoContainer *hcal_towers = findNode::getClass<TowerInfoContainer>(topNode, IHTowerName);
  if (!hcal_towers)
  {
    std::cout << PHWHERE << "Cannot find " << IHTowerName << " for creating new tower containers. Exiting" << std::endl;
    exit(1);
  }

  // store the new EMCal towers

  PHCompositeNode *emcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "CEMC"));
  if (!emcalNode)
  {
    std::cout << PHWHERE << "EMCal Node note found, doing nothing." << std::endl;
  }
  EMTowerName = emcal_output_node();
  TowerInfoContainer *test_emcal_tower = findNode::getClass<TowerInfoContainer>(topNode, EMTowerName);
  if (!test_emcal_tower)
  {
    if (Verbosity() > 0)
    {
      std::cout << "SubtractTowers::CreateNode : creating " << EMTowerName << " node " << std::endl;
    }

    TowerInfoContainer *emcal_towers = dynamic_cast<TowerInfoContainer *>(hcal_towers->CloneMe());
    PHIODataNode<PHObject> *emcalTowerNode = new PHIODataNode<PHObject>(emcal_towers, EMTowerName, "PHObject");
    emcalNode->addNode(emcalTowerNode);
  }
  else
  {
    std::cout << "SubtractTowers::CreateNode : " << EMTowerName << " already exists! " << std::endl;
  }

  // store the new IHCal towers
  PHCompositeNode *ihcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "HCALIN"));
  if (!ihcalNode)
  {
    std::cout << PHWHERE << "IHCal Node note found, doing nothing." << std::endl;
  }
  IHTowerName = ihcal_output_node();
  TowerInfoContainer *test_ihcal_tower = findNode::getClass<TowerInfoContainer>(topNode, IHTowerName);
  if (!test_ihcal_tower)
  {
    if (Verbosity() > 0)
    {
      std::cout << "SubtractTowers::CreateNode : creating " << IHTowerName << " node " << std::endl;
    }

    TowerInfoContainer *ihcal_towers = dynamic_cast<TowerInfoContainer *>(hcal_towers->CloneMe());
    PHIODataNode<PHObject> *ihcalTowerNode = new PHIODataNode<PHObject>(ihcal_towers, IHTowerName, "PHObject");
    ihcalNode->addNode(ihcalTowerNode);
  }
  else
  {
    std::cout << "SubtractTowers::CreateNode : " << IHTowerName << " already exists! " << std::endl;
  }

  // store the new OHCal towers
  PHCompositeNode *ohcalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "HCALOUT"));
  if (!ohcalNode)
  {
    std::cout << PHWHERE << "OHCal Node note found, doing nothing." << std::endl;
  }
  OHTowerName = ohcal_output_node();
  TowerInfoContainer *test_ohcal_tower = findNode::getClass<TowerInfoContainer>(topNode, OHTowerName);
  if (!test_ohcal_tower)
  {
    if (Verbosity() > 0)
    {
      std::cout << "SubtractTowers::CreateNode : creating " << OHTowerName << " node " << std::endl;
    }

    TowerInfoContainer *ohcal_towers = dynamic_cast<TowerInfoContainer *>(hcal_towers->CloneMe());
    PHIODataNode<PHObject> *ohcalTowerNode = new PHIODataNode<PHObject>(ohcal_towers, OHTowerName, "PHObject");
    ohcalNode->addNode(ohcalTowerNode);
  }
  else
  {
    std::cout << "SubtractTowers::CreateNode : " << OHTowerName << " already exists! " << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
