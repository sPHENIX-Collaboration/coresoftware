#include "TpcConditionsReco.h"
#include "TpcConditions.h"

#include <ffamodules/CDBInterface.h>

#include <cdbobjects/CDBTTree.h>
#include <algorithm>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <ffarawobjects/Gl1Packet.h>

#include <cstdint>
#include <iostream>

#include <string>
#include <vector>

//  This namespace organizes Charles' "Undamaged GEM List"
//  into vectors for simpler consumption later.
namespace
{
  const std::vector<std::string> SR1 = {
    "S_01_R1_G4_IMon",
    "S_04_R1_G4_IMon",
    "S_05_R1_G4_IMon",
    "S_10_R1_G4_IMon"
  };
  
  const std::vector<std::string> SR2 = {
    "S_01_R2_G4_IMon",
    "S_02_R2_G4_IMon",
    "S_07_R2_G4_IMon",
    "S_09_R2_G4_IMon",
    "S_10_R2_G4_IMon",
    "S_12_R2_G4_IMon"
  };
  
  const std::vector<std::string> SR3 = {
    "S_08_R3_G4_IMon",
    "S_09_R3_G4_IMon",
    "S_10_R3_G4_IMon",
  };
  
  const std::vector<std::string> NR1 = {
    "N_03_R1_G4_IMon",
    "N_04_R1_G4_IMon",
    "N_06_R1_G4_IMon",
    "N_07_R1_G4_IMon",
    "N_08_R1_G4_IMon",
    "N_10_R1_G4_IMon",
    "N_11_R1_G4_IMon",
    "N_12_R1_G4_IMon"
  };
  
  const std::vector<std::string> NR2 = {
    "N_02_R2_G4_IMon",
    "N_04_R2_G4_IMon",
    "N_05_R2_G4_IMon",
    "N_09_R2_G4_IMon",
    "N_10_R2_G4_IMon",
    "N_11_R2_G4_IMon",
    "N_12_R2_G4_IMon"
  };
  
  const std::vector<std::string> NR3 = {
    "N_02_R3_G4_IMon",
    "N_10_R3_G4_IMon"
  };
  
  const std::vector<std::string> SOUTH = []()
  {
    std::vector<std::string> v;
    v.insert(v.end(), SR1.begin(), SR1.end());
    v.insert(v.end(), SR2.begin(), SR2.end());
    v.insert(v.end(), SR3.begin(), SR3.end());
    return v;
  }();
  
  const std::vector<std::string> NORTH = []()
  {
    std::vector<std::string> v;
    v.insert(v.end(), NR1.begin(), NR1.end());
    v.insert(v.end(), NR2.begin(), NR2.end());
    v.insert(v.end(), NR3.begin(), NR3.end());
    return v;
  }();
  
  const std::vector<std::string> ALL = []()
  {
    std::vector<std::string> v;
    v.insert(v.end(), SOUTH.begin(), SOUTH.end());
    v.insert(v.end(), NORTH.begin(), NORTH.end());
    return v;
  }();
}

TpcConditionsReco::TpcConditionsReco(const std::string &name)
  : SubsysReco(name)
{
}

TpcConditionsReco::~TpcConditionsReco()
{
  delete m_tree;
}

int TpcConditionsReco::InitRun(PHCompositeNode *topNode)
{
  std::cout << "TpcConditionsReco::InitRun(PHCompositeNode *topNode) Initializing"
            << std::endl;

  // Get the TPC conditions payload from CDB
  m_cdb = CDBInterface::instance();
  std::string calibdir = m_cdb->getUrl("TPC_CONDITIONS");
  m_tree = new CDBTTree(calibdir);
  m_tree->LoadCalibrations();

  // Build the BCO -> CDBTTree channel lookup
  m_bco_to_channel.clear();
  for (unsigned int channel = 0;
       channel < m_tree->GetUInt64EntryMap().size();
       ++channel)
  {
    uint64_t bco = m_tree->GetUInt64Value(channel, "bco");
    m_bco_to_channel[bco] = channel;
  }

  // Find/create RUN node
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode =
      dynamic_cast<PHCompositeNode *>(
          iter.findFirst("PHCompositeNode", "RUN"));

  if (!runNode)
  {
    runNode = new PHCompositeNode("RUN");
    topNode->addNode(runNode);
  }

  // Create the transient TPC conditions object
  m_conditions = new TpcConditions();
  PHDataNode<PHObject> *conditionsNode =
      new PHDataNode<PHObject>(
          m_conditions, "TpcConditions", "PHObject");
  runNode->addNode(conditionsNode);

  return Fun4AllReturnCodes::EVENT_OK;
}


int TpcConditionsReco::process_event(PHCompositeNode *topNode)
{
  // Get the current event BCO
  Gl1Packet *gl1 =
      findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1)
  {
    std::cout << "TpcConditionsReco::process_event - "
              << "could not find GL1RAWHIT node"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  uint64_t bco = gl1->getBCO();

  // Find the most recent conditions sample at or before this BCO
  auto iter = m_bco_to_channel.upper_bound(bco);
  if (iter == m_bco_to_channel.begin())
  {
    std::cout << "TpcConditionsReco::process_event - "
              << "no TPC conditions available for BCO "
              << bco << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  --iter;
  int channel = iter->second;
  
  // Fill current-event conditions
  m_conditions->set_Temperature(m_tree->GetFloatValue(channel, "gas_temperature"));
  m_conditions->set_Pressure   (m_tree->GetFloatValue(channel, "gas_pressure"));
  m_conditions->set_FieldOK    (m_tree->GetFloatValue(channel, "FieldOK") != 0.0);
  m_conditions->set_GainOK     (m_tree->GetFloatValue(channel, "GainOK") != 0.0);

  // Here are the critical LoadCurrent calculations.
  m_conditions->set_LoadCurrent(get_MedianCurrent(channel, ALL));
  m_conditions->set_LoadNorth  (get_MedianCurrent(channel, NORTH));
  m_conditions->set_LoadSouth  (get_MedianCurrent(channel, SOUTH));
  m_conditions->set_LoadSR1    (get_MedianCurrent(channel, SR1));
  m_conditions->set_LoadSR2    (get_MedianCurrent(channel, SR2));
  m_conditions->set_LoadSR3    (get_MedianCurrent(channel, SR3));
  m_conditions->set_LoadNR1    (get_MedianCurrent(channel, NR1));
  m_conditions->set_LoadNR2    (get_MedianCurrent(channel, NR2));
  m_conditions->set_LoadNR3    (get_MedianCurrent(channel, NR3));

  // Do or die
  if (!m_conditions->get_FieldOK() ||
      !m_conditions->get_GainOK())
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

float TpcConditionsReco::get_MedianCurrent(
    int channel,
    const std::vector<std::string> &channels)
{
  std::vector<float> currents;
  currents.reserve(channels.size());

  for (const auto &name : channels)
  {
    currents.push_back(
        m_tree->GetFloatValue(channel, name));
  }

  std::sort(currents.begin(), currents.end());

  const size_t n = currents.size();

  if (n % 2)
  {
    return currents[n / 2];
  }

  return 0.5F *
         (currents[n / 2 - 1] + currents[n / 2]);
}
