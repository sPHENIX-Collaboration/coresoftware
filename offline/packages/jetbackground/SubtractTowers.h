#ifndef JETBACKGROUND_SUBTRACTTOWERS_H
#define JETBACKGROUND_SUBTRACTTOWERS_H

//===========================================================
/// \file SubtractTowers.h
/// \brief creates new UE-subtracted towers
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>

// forward declarations
class PHCompositeNode;

/// \class SubtractTowers
///
/// \brief creates new UE-subtracted towers
///
/// Using a previously determined background UE density, this module
/// constructs a new set of towers by subtracting the background from
/// existing raw towers
///
class SubtractTowers : public SubsysReco
{
 public:
  SubtractTowers(const std::string &name = "SubtractTowers");
  ~SubtractTowers() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void SetFlowModulation(bool use_flow_modulation) { _use_flow_modulation = use_flow_modulation; }
  void set_towerinfo(bool /*use_towerinfo*/)
  {
    // m_use_towerinfo = use_towerinfo;
  }
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
  // output tower nodes; empty (default) means
  // m_towerNodePrefix + "_CEMC_RETOWER_SUB1" / "_HCALIN_SUB1" / "_HCALOUT_SUB1"
  void set_emcal_output_node(const std::string &name) { m_emcal_output_node = name; }
  void set_ihcal_output_node(const std::string &name) { m_ihcal_output_node = name; }
  void set_ohcal_output_node(const std::string &name) { m_ohcal_output_node = name; }
  void set_inputTowerBackgroundNode(const std::string &nodeName)
  {
    m_towerBackgroundNode = nodeName;
  }

 private:
  int CreateNode(PHCompositeNode *topNode);

  bool _use_flow_modulation{false};
  std::string m_towerNodePrefix{"TOWERINFO_CALIB"};
  std::string m_towerBackgroundNode{"TowerInfoBackground_Sub2"};
  std::string m_emcal_input_node;
  std::string m_ihcal_input_node;
  std::string m_ohcal_input_node;
  std::string m_ihcal_geom_node{"TOWERGEOM_HCALIN"};
  std::string m_ohcal_geom_node{"TOWERGEOM_HCALOUT"};
  std::string emcal_input_node() const { return m_emcal_input_node.empty() ? m_towerNodePrefix + "_CEMC_RETOWER" : m_emcal_input_node; }
  std::string ihcal_input_node() const { return m_ihcal_input_node.empty() ? m_towerNodePrefix + "_HCALIN" : m_ihcal_input_node; }
  std::string ohcal_input_node() const { return m_ohcal_input_node.empty() ? m_towerNodePrefix + "_HCALOUT" : m_ohcal_input_node; }
  std::string m_emcal_output_node;
  std::string m_ihcal_output_node;
  std::string m_ohcal_output_node;
  std::string emcal_output_node() const { return m_emcal_output_node.empty() ? m_towerNodePrefix + "_CEMC_RETOWER_SUB1" : m_emcal_output_node; }
  std::string ihcal_output_node() const { return m_ihcal_output_node.empty() ? m_towerNodePrefix + "_HCALIN_SUB1" : m_ihcal_output_node; }
  std::string ohcal_output_node() const { return m_ohcal_output_node.empty() ? m_towerNodePrefix + "_HCALOUT_SUB1" : m_ohcal_output_node; }
  std::string EMTowerName;
  std::string IHTowerName;
  std::string OHTowerName;
};

#endif
