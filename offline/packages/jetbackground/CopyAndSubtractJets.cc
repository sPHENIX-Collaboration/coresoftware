#include "CopyAndSubtractJets.h"

#include "TowerBackground.h"

#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>

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
#include <utility>
#include <vector>

CopyAndSubtractJets::CopyAndSubtractJets(const std::string &name)
  : SubsysReco(name)
{
}

int CopyAndSubtractJets::InitRun(PHCompositeNode *topNode)
{
  CreateNode(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CopyAndSubtractJets::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "CopyAndSubtractJets::process_event: entering, with _use_flow_modulation = " << _use_flow_modulation << std::endl;
  }

  // pull out needed calo tower info
  EMTowerName = emcal_input_node();
  IHTowerName = ihcal_input_node();
  OHTowerName = ohcal_input_node();
  TowerInfoContainer *towerinfosEM3 = findNode::getClass<TowerInfoContainer>(topNode, EMTowerName);
  TowerInfoContainer *towerinfosIH3 = findNode::getClass<TowerInfoContainer>(topNode, IHTowerName);
  TowerInfoContainer *towerinfosOH3 = findNode::getClass<TowerInfoContainer>(topNode, OHTowerName);
  if (!towerinfosEM3)
  {
    std::cout << "CopyAndSubtractJets::process_event: Cannot find node " << EMTowerName << std::endl;
    exit(1);
  }
  if (!towerinfosIH3)
  {
    std::cout << "CopyAndSubtractJets::process_event: Cannot find node " << IHTowerName << std::endl;
    exit(1);
  }
  if (!towerinfosOH3)
  {
    std::cout << "CopyAndSubtractJets::process_event: Cannot find node " << OHTowerName << std::endl;
    exit(1);
  }

  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, m_ihcal_geom_node);
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, m_ohcal_geom_node);

  // pull out jets and background
  JetContainer *unsub_jets = findNode::getClass<JetContainer>(topNode, m_rawseed_node);
  JetContainer *sub_jets = findNode::getClass<JetContainer>(topNode, m_subseed_node);
  TowerBackground *background = findNode::getClass<TowerBackground>(topNode, m_background_node);
  std::vector<float> background_UE_0 = background->get_UE(0);
  std::vector<float> background_UE_1 = background->get_UE(1);
  std::vector<float> background_UE_2 = background->get_UE(2);

  float background_v2 = background->get_v2();
  float background_Psi2 = background->get_Psi2();

  if (Verbosity() > 0)
  {
    std::cout << "CopyAndSubtractJets::process_event: entering with # unsubtracted jets = " << unsub_jets->size() << std::endl;
    std::cout << "CopyAndSubtractJets::process_event: entering with # subtracted jets = " << sub_jets->size() << std::endl;
  }

  // iterate over old jets
  int ijet = 0;
  for (auto *this_jet : *unsub_jets)
  {
    float this_pt = this_jet->get_pt();
    float this_phi = this_jet->get_phi();
    float this_eta = this_jet->get_eta();

    float new_total_px = 0;
    float new_total_py = 0;
    float new_total_pz = 0;
    float new_total_e = 0;

    // if (this_jet->get_pt() < 5) continue;

    if (Verbosity() > 1 && this_jet->get_pt() > 5)
    {
      std::cout << "CopyAndSubtractJets::process_event: unsubtracted jet with pt / eta / phi = " << this_pt << " / " << this_eta << " / " << this_phi << std::endl;
    }

    for (const auto &comp : this_jet->get_comp_vec())
    {
      RawTowerGeom *tower_geom = nullptr;
      TowerInfo *towerinfo = nullptr;

      double comp_e = 0;
      double comp_eta = 0;
      double comp_phi = 0;

      int comp_ieta = 0;

      double comp_background = 0;

      if (comp.first == Jet::SRC::HCALIN_TOWER || comp.first == Jet::SRC::HCALIN_TOWERINFO)
      {
        towerinfo = towerinfosIH3->get_tower_at_channel(comp.second);
        unsigned int towerkey = towerinfosIH3->encode_key(comp.second);
        comp_ieta = towerinfosIH3->getTowerEtaBin(towerkey);
        int comp_iphi = towerinfosIH3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, comp_ieta, comp_iphi);

        tower_geom = geomIH->get_tower_geometry(key);
        comp_background = background_UE_1.at(comp_ieta);
      }
      else if (comp.first == Jet::SRC::HCALOUT_TOWER || comp.first == Jet::SRC::HCALOUT_TOWERINFO)
      {
        towerinfo = towerinfosOH3->get_tower_at_channel(comp.second);
        unsigned int towerkey = towerinfosOH3->encode_key(comp.second);
        comp_ieta = towerinfosOH3->getTowerEtaBin(towerkey);
        int comp_iphi = towerinfosOH3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, comp_ieta, comp_iphi);
        tower_geom = geomOH->get_tower_geometry(key);
        comp_background = background_UE_2.at(comp_ieta);
      }
      else if (comp.first == Jet::SRC::CEMC_TOWER_RETOWER || comp.first == Jet::SRC::CEMC_TOWERINFO_RETOWER)
      {
        towerinfo = towerinfosEM3->get_tower_at_channel(comp.second);
        unsigned int towerkey = towerinfosEM3->encode_key(comp.second);
        comp_ieta = towerinfosEM3->getTowerEtaBin(towerkey);
        int comp_iphi = towerinfosEM3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, comp_ieta, comp_iphi);

        tower_geom = geomIH->get_tower_geometry(key);
        comp_background = background_UE_0.at(comp_ieta);
      }
      if (towerinfo)
      {
        comp_e = towerinfo->get_energy();
      }

      if (tower_geom)
      {
        comp_eta = tower_geom->get_eta();
        comp_phi = tower_geom->get_phi();
      }

      if (Verbosity() > 4 && this_jet->get_pt() > 5)
      {
        std::cout << "CopyAndSubtractJets::process_event: --> constituent in layer " << comp.first << ", has unsub E = " << comp_e << ", is at ieta #" << comp_ieta << ", and has UE = " << comp_background << std::endl;
      }

      // flow modulate background if turned on
      if (_use_flow_modulation)
      {
        comp_background = comp_background * (1 + 2 * background_v2 * cos(2 * (comp_phi - background_Psi2)));
        if (Verbosity() > 4 && this_jet->get_pt() > 5)
        {
          std::cout << "CopyAndSubtractJets::process_event: --> --> flow mod, at phi = " << comp_phi << ", v2 and Psi2 are = " << background_v2 << " , " << background_Psi2 << ", UE after modulation = " << comp_background << std::endl;
        }
      }

      // update constituent energy based on the background
      double comp_sub_e = comp_e - comp_background;

      // now define new kinematics

      double comp_px = comp_sub_e / cosh(comp_eta) * cos(comp_phi);
      double comp_py = comp_sub_e / cosh(comp_eta) * sin(comp_phi);
      double comp_pz = comp_sub_e * tanh(comp_eta);

      new_total_px += comp_px;
      new_total_py += comp_py;
      new_total_pz += comp_pz;
      new_total_e += comp_sub_e;
    }

    auto *new_jet = sub_jets->add_jet();  // returns a new Jet_v2

    new_jet->set_px(new_total_px);
    new_jet->set_py(new_total_py);
    new_jet->set_pz(new_total_pz);
    new_jet->set_e(new_total_e);
    new_jet->set_id(ijet);

    if (Verbosity() > 1 && this_pt > 5)
    {
      std::cout << "CopyAndSubtractJets::process_event: old jet #" << ijet << ", old px / py / pz / e = " << this_jet->get_px() << " / " << this_jet->get_py() << " / " << this_jet->get_pz() << " / " << this_jet->get_e() << std::endl;
      std::cout << "CopyAndSubtractJets::process_event: new jet #" << ijet << ", new px / py / pz / e = " << new_jet->get_px() << " / " << new_jet->get_py() << " / " << new_jet->get_pz() << " / " << new_jet->get_e() << std::endl;
    }

    ijet++;
  }

  if (Verbosity() > 0)
  {
    std::cout << "CopyAndSubtractJets::process_event: exiting with # subtracted jets = " << sub_jets->size() << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CopyAndSubtractJets::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int CopyAndSubtractJets::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Looking for the ANTIKT node
  PHCompositeNode *antiktNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", m_jet_node));
  if (!antiktNode)
  {
    std::cout << PHWHERE << m_jet_node << " node not found, doing nothing." << std::endl;
  }

  // Looking for the TOWER node
  PHCompositeNode *towerNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", m_input_node));
  if (!towerNode)
  {
    std::cout << PHWHERE << m_input_node << " node not found, doing nothing." << std::endl;
  }

  // store the new jet collection
  JetContainer *test_jets = findNode::getClass<JetContainer>(topNode, m_subseed_node);
  if (!test_jets)
  {
    if (Verbosity() > 0)
    {
      std::cout << "CopyAndSubtractJets::CreateNode : creating " << m_subseed_node << " node " << std::endl;
    }
    JetContainer *sub_jets = new JetContainerv1();
    PHIODataNode<PHObject> *subjetNode = new PHIODataNode<PHObject>(sub_jets, m_subseed_node, "PHObject");
    towerNode->addNode(subjetNode);
  }
  else
  {
    std::cout << "CopyAndSubtractJets::CreateNode : " << m_subseed_node << " already exists! " << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
