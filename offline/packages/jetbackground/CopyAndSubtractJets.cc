#include "CopyAndSubtractJets.h"

#include "TowerBackground.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
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
  RawTowerContainer *towersEM3 = nullptr;
  RawTowerContainer *towersIH3 = nullptr;
  RawTowerContainer *towersOH3 = nullptr;
  TowerInfoContainer *towerinfosEM3 = nullptr;
  TowerInfoContainer *towerinfosIH3 = nullptr;
  TowerInfoContainer *towerinfosOH3 = nullptr;
  if (m_use_towerinfo)
  {
    towerinfosEM3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
    towerinfosIH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
    towerinfosOH3 = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  }
  else
  {
    towersEM3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC_RETOWER");
    towersIH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
    towersOH3 = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");
  }

  RawTowerGeomContainer *geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  RawTowerGeomContainer *geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  // pull out jets and background
  JetContainer *unsub_jets;
  JetContainer *sub_jets;
  TowerBackground *background;
  if (m_use_towerinfo)
  {
    unsub_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
    sub_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_HIRecoSeedsSub_r02");
    background = findNode::getClass<TowerBackground>(topNode, "TowerInfoBackground_Sub1");
  }
  else
  {
    unsub_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Tower_HIRecoSeedsRaw_r02");
    sub_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Tower_HIRecoSeedsSub_r02");
    background = findNode::getClass<TowerBackground>(topNode, "TowerBackground_Sub1");
  }
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
  for (auto this_jet : *unsub_jets)
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
      RawTower *tower = nullptr;
      RawTowerGeom *tower_geom = nullptr;
      TowerInfo *towerinfo = nullptr;

      double comp_e = 0;
      double comp_eta = 0;
      double comp_phi = 0;

      int comp_ieta = 0;

      double comp_background = 0;

      if (m_use_towerinfo)
      {
        if (comp.first == 5 || comp.first == 26)
        {
          towerinfo = towerinfosIH3->get_tower_at_channel(comp.second);
          unsigned int towerkey = towerinfosIH3->encode_key(comp.second);
          comp_ieta = towerinfosIH3->getTowerEtaBin(towerkey);
          int comp_iphi = towerinfosIH3->getTowerPhiBin(towerkey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, comp_ieta, comp_iphi);

          tower_geom = geomIH->get_tower_geometry(key);
          comp_background = background_UE_1.at(comp_ieta);
        }
        else if (comp.first == 7 || comp.first == 27)
        {
          towerinfo = towerinfosOH3->get_tower_at_channel(comp.second);
          unsigned int towerkey = towerinfosOH3->encode_key(comp.second);
          comp_ieta = towerinfosOH3->getTowerEtaBin(towerkey);
          int comp_iphi = towerinfosOH3->getTowerPhiBin(towerkey);
          const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, comp_ieta, comp_iphi);
          tower_geom = geomOH->get_tower_geometry(key);
          comp_background = background_UE_2.at(comp_ieta);
        }
        else if (comp.first == 13 || comp.first == 28)
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
      }
      else
      {
        if (comp.first == 5)
        {
          tower = towersIH3->getTower(comp.second);
          tower_geom = geomIH->get_tower_geometry(tower->get_key());

          comp_ieta = geomIH->get_etabin(tower_geom->get_eta());
          comp_background = background_UE_1.at(comp_ieta);
        }
        else if (comp.first == 7)
        {
          tower = towersOH3->getTower(comp.second);
          tower_geom = geomOH->get_tower_geometry(tower->get_key());

          comp_ieta = geomOH->get_etabin(tower_geom->get_eta());
          comp_background = background_UE_2.at(comp_ieta);
        }
        else if (comp.first == 13)
        {
          tower = towersEM3->getTower(comp.second);
          tower_geom = geomIH->get_tower_geometry(tower->get_key());

          comp_ieta = geomIH->get_etabin(tower_geom->get_eta());
          comp_background = background_UE_0.at(comp_ieta);
        }
        if (tower)
        {
          comp_e = tower->get_energy();
        }
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

    auto new_jet = sub_jets->add_jet();  // returns a new Jet_v2

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
  PHCompositeNode *antiktNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "ANTIKT"));
  if (!antiktNode)
  {
    std::cout << PHWHERE << "ANTIKT node not found, doing nothing." << std::endl;
  }

  // Looking for the TOWER node
  PHCompositeNode *towerNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "TOWER"));
  if (!towerNode)
  {
    std::cout << PHWHERE << "TOWER node not found, doing nothing." << std::endl;
  }

  // store the new jet collection
  JetContainer *test_jets;
  if (m_use_towerinfo)
  {
    test_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_HIRecoSeedsSub_r02");
  }
  else
  {
    test_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Tower_HIRecoSeedsSub_r02");
  }
  if (!test_jets)
  {
    JetContainer *sub_jets = new JetContainerv1();
    PHIODataNode<PHObject> *subjetNode;
    if (m_use_towerinfo)
    {
      if (Verbosity() > 0)
      {
        std::cout << "CopyAndSubtractJets::CreateNode : creating AntiKt_TowerInfo_HIRecoSeedsSub_r02 node " << std::endl;
      }
      subjetNode = new PHIODataNode<PHObject>(sub_jets, "AntiKt_TowerInfo_HIRecoSeedsSub_r02", "PHObject");
    }
    else
    {
      if (Verbosity() > 0)
      {
        std::cout << "CopyAndSubtractJets::CreateNode : creating AntiKt_Tower_HIRecoSeedsSub_r02 node " << std::endl;
      }
      subjetNode = new PHIODataNode<PHObject>(sub_jets, "AntiKt_Tower_HIRecoSeedsSub_r02", "PHObject");
    }
    towerNode->addNode(subjetNode);
  }
  else
  {
    std::cout << "CopyAndSubtractJets::CreateNode : AntiKt_Tower_HIRecoSeedsSub_r02 already exists! " << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
