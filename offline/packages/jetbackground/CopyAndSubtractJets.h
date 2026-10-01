#ifndef JETBACKGROUND_COPYANDSUBTRACTJETS_H
#define JETBACKGROUND_COPYANDSUBTRACTJETS_H

//===========================================================
/// \file CopyAndSubtractJets.h
/// \brief Creates subtracted copy of a jet collection
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>

// forward declarations
class PHCompositeNode;

/// \class CopyAndSubtractJets
///
/// \brief Creates subtractd copy of a jet collection
///
/// Makes a copy of a jet collection with a new name and then updates
/// the kinematics of the jets in that collection based on a given UE
/// background (intended use is to create the set of jets used as
/// seeds in the second part of UE determination procedure)
///
class CopyAndSubtractJets : public SubsysReco
{
 public:
  CopyAndSubtractJets(const std::string &name = "CopyAndSubtractJets");
  ~CopyAndSubtractJets() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void SetFlowModulation(bool use_flow_modulation) { _use_flow_modulation = use_flow_modulation; }
  // kept so existing macros still compile; only TowerInfo input is supported
  void set_towerinfo(bool /*use_towerinfo*/) {}
  void set_towerNodePrefix(const std::string &prefix)
  {
    m_towerNodePrefix = prefix;
    return;
  }
  // input tower / geometry nodes; an empty tower node (default) means
  // m_towerNodePrefix + "_CEMC_RETOWER" / "_HCALIN" / "_HCALOUT"
  void set_emcal_input_node(const std::string &name) { m_emcal_input_node = name; }
  void set_ihcal_input_node(const std::string &name) { m_ihcal_input_node = name; }
  void set_ohcal_input_node(const std::string &name) { m_ohcal_input_node = name; }
  void set_ihcal_geom_node(const std::string &name) { m_ihcal_geom_node = name; }
  void set_ohcal_geom_node(const std::string &name) { m_ohcal_geom_node = name; }

  void set_rawseed_node(const std::string &name) { m_rawseed_node = name; }
  void set_subseed_node(const std::string &name) { m_subseed_node = name; }
  void set_background_node(const std::string &name) { m_background_node = name; }
  void set_jet_node(const std::string &name) { m_jet_node = name; }
  void set_input_node(const std::string &name) { m_input_node = name; }

 private:
  int CreateNode(PHCompositeNode *topNode);

  bool _use_flow_modulation{false};
  std::string m_towerNodePrefix{"TOWERINFO_CALIB"};
  std::string m_emcal_input_node{};
  std::string m_ihcal_input_node{};
  std::string m_ohcal_input_node{};
  std::string m_ihcal_geom_node{"TOWERGEOM_HCALIN"};
  std::string m_ohcal_geom_node{"TOWERGEOM_HCALOUT"};
  std::string emcal_input_node() const { return m_emcal_input_node.empty() ? m_towerNodePrefix + "_CEMC_RETOWER" : m_emcal_input_node; }
  std::string ihcal_input_node() const { return m_ihcal_input_node.empty() ? m_towerNodePrefix + "_HCALIN" : m_ihcal_input_node; }
  std::string ohcal_input_node() const { return m_ohcal_input_node.empty() ? m_towerNodePrefix + "_HCALOUT" : m_ohcal_input_node; }
  std::string EMTowerName;
  std::string IHTowerName;
  std::string OHTowerName;

  std::string m_rawseed_node = "AntiKt_TowerInfo_HIRecoSeedsRaw_r02";
  std::string m_subseed_node = "AntiKt_TowerInfo_HIRecoSeedsSub_r02";
  std::string m_background_node = "TowerInfoBackground_Sub1";
  std::string m_jet_node = "ANTIKT";
  std::string m_input_node = "TOWER";
};

#endif
