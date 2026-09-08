#include "TpcConditionsReco.h"
#include "TpcConditions.h"

#include <ffamodules/CDBInterface.h>

#include <cdbobjects/CDBTTree.h>

#include <ffarawobjects/Gl1Packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <phool/recoConsts.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
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
      "S_10_R1_G4_IMon"};

  const std::vector<std::string> SR2 = {
      "S_01_R2_G4_IMon",
      "S_02_R2_G4_IMon",
      "S_07_R2_G4_IMon",
      "S_09_R2_G4_IMon",
      "S_10_R2_G4_IMon",
      "S_12_R2_G4_IMon"};

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
      "N_12_R1_G4_IMon"};

  const std::vector<std::string> NR2 = {
      "N_02_R2_G4_IMon",
      "N_04_R2_G4_IMon",
      "N_05_R2_G4_IMon",
      "N_09_R2_G4_IMon",
      "N_10_R2_G4_IMon",
      "N_11_R2_G4_IMon",
      "N_12_R2_G4_IMon"};

  const std::vector<std::string> NR3 = {
      "N_02_R3_G4_IMon",
      "N_10_R3_G4_IMon"};

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
}  // namespace

TpcConditionsReco::TpcConditionsReco(const std::string &name)
  : SubsysReco(name)
{
}

TpcConditionsReco::~TpcConditionsReco()
{
  delete m_tree;
}

void TpcConditionsReco::fillConditions(int channel)
{
  m_conditions->set_Temperature(m_tree->GetFloatValue(channel, "gas_temperature"));
  m_conditions->set_Pressure(m_tree->GetFloatValue(channel, "gas_pressure"));
  m_conditions->set_FieldOK(m_tree->GetFloatValue(channel, "FieldOK") != 0.0);
  m_conditions->set_GainOK(m_tree->GetFloatValue(channel, "GainOK") != 0.0);

  m_conditions->set_LoadCurrent(get_MedianCurrent(channel, ALL));
  m_conditions->set_LoadNorth(get_MedianCurrent(channel, NORTH));
  m_conditions->set_LoadSouth(get_MedianCurrent(channel, SOUTH));
  m_conditions->set_LoadSR1(get_MedianCurrent(channel, SR1));
  m_conditions->set_LoadSR2(get_MedianCurrent(channel, SR2));
  m_conditions->set_LoadSR3(get_MedianCurrent(channel, SR3));
  m_conditions->set_LoadNR1(get_MedianCurrent(channel, NR1));
  m_conditions->set_LoadNR2(get_MedianCurrent(channel, NR2));
  m_conditions->set_LoadNR3(get_MedianCurrent(channel, NR3));
}

float TpcConditionsReco::get_AverageMedianCurrent(
    const std::vector<std::string> &channels)
{
  if (m_bco_to_channel.empty())
  {
    return 0.0F;
  }

  double sum = 0.0;

  for (const auto &[bco, channel] : m_bco_to_channel)
  {
    sum += get_MedianCurrent(channel, channels);
  }

  return static_cast<float>(sum / static_cast<double>(m_bco_to_channel.size()));
}

int TpcConditionsReco::InitRun(PHCompositeNode *topNode)
{
  std::cout << "TpcConditionsReco::InitRun(PHCompositeNode *topNode) Initializing"
            << std::endl;

  Gl1Packet *gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1)
  {
    gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
  }

  if (!gl1)
  {
    std::cout << "TpcConditionsReco::InitRun - could not find GL1 packet"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const int64_t rawBco = gl1->lValue(0, "BCO");

  if (rawBco <= 0)
  {
    std::cout << "TpcConditionsReco::InitRun - invalid GL1 BCO "
              << rawBco << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  const uint64_t targetBco = static_cast<uint64_t>(rawBco);

  std::cout << "\033[1;31m"
            << " targetBCO = "
            << targetBco
            << "\033[0m" << std::endl;

  // Get the TPC conditions payload from CDB
  std::string calibdir = CDBInterface::instance()->getUrl("TPC_CONDITIONS");
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
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));

  if (!parNode)
  {
    parNode = new PHCompositeNode("PAR");
    topNode->addNode(parNode);
  }

  // Create the transient TPC conditions object
  m_conditions = new TpcConditions();
  PHDataNode<TpcConditions> *conditionsNode = new PHDataNode<TpcConditions>(m_conditions, "TpcConditions", "Data");
  parNode->addNode(conditionsNode);

  m_conditions->set_AverageLoadCurrent(get_AverageMedianCurrent(ALL));
  m_conditions->set_AverageLoadNorth(get_AverageMedianCurrent(NORTH));
  m_conditions->set_AverageLoadSouth(get_AverageMedianCurrent(SOUTH));

  m_conditions->set_AverageLoadSR1(get_AverageMedianCurrent(SR1));
  m_conditions->set_AverageLoadSR2(get_AverageMedianCurrent(SR2));
  m_conditions->set_AverageLoadSR3(get_AverageMedianCurrent(SR3));

  m_conditions->set_AverageLoadNR1(get_AverageMedianCurrent(NR1));
  m_conditions->set_AverageLoadNR2(get_AverageMedianCurrent(NR2));
  m_conditions->set_AverageLoadNR3(get_AverageMedianCurrent(NR3));

  if (!m_bco_to_channel.empty())
  {
    recoConsts *rc = recoConsts::instance();
    const int segment = rc->get_IntFlag("RUNSEGMENT");

    // Find closest stored conditions BCO.
    auto upper = m_bco_to_channel.upper_bound(targetBco);
    auto selected = upper;

    if (upper == m_bco_to_channel.end())
    {
      selected = std::prev(m_bco_to_channel.end());
    }
    else if (upper != m_bco_to_channel.begin())
    {
      auto lower = std::prev(upper);

      if ((targetBco - lower->first) <= (upper->first - targetBco))
      {
        selected = lower;
      }
    }

    if (!m_bco_to_channel.empty())
    {
      const auto first = m_bco_to_channel.begin();
      const auto last = std::prev(m_bco_to_channel.end());
      const auto middle = std::next(
          m_bco_to_channel.begin(),
          static_cast<long>(m_bco_to_channel.size() / 2));

      std::cout << "\033[1;31m"
                << "R1 conditions:"
                << "\n  FIRST  BCO=" << first->first
                << " NR1=" << get_MedianCurrent(first->second, NR1)
                << " SR1=" << get_MedianCurrent(first->second, SR1)
                << "\n  MIDDLE BCO=" << middle->first
                << " NR1=" << get_MedianCurrent(middle->second, NR1)
                << " SR1=" << get_MedianCurrent(middle->second, SR1)
                << "\n  LAST   BCO=" << last->first
                << " NR1=" << get_MedianCurrent(last->second, NR1)
                << " SR1=" << get_MedianCurrent(last->second, SR1)
                << "\n  AVERAGE NR1=" << m_conditions->get_AverageLoadNR1()
                << " SR1=" << m_conditions->get_AverageLoadSR1()
                << "\033[0m" << std::endl;
    }

    fillConditions(selected->second);

    std::cout << "TpcConditionsReco::InitRun - segment " << segment
              << ", target BCO " << targetBco
              << ", conditions BCO " << selected->first
              << ", LoadNR1 " << m_conditions->get_LoadNR1()
              << ", LoadSR1 " << m_conditions->get_LoadSR1()
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcConditionsReco::process_event(PHCompositeNode *topNode)
{
  // Get the current event BCO
  Gl1Packet *gl1 = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1)
  {
    std::cout << "TpcConditionsReco::process_event - "
              << "could not find GL1RAWHIT node"
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  uint64_t bco = gl1->getBCO();

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

  m_conditions->set_FieldOK(m_tree->GetFloatValue(channel, "FieldOK") != 0.0);
  m_conditions->set_GainOK(m_tree->GetFloatValue(channel, "GainOK") != 0.0);


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
    currents.push_back(m_tree->GetFloatValue(channel, name));
  }

  std::sort(currents.begin(), currents.end());

  const size_t n = currents.size();

  if (n % 2)
  {
    return currents[n / 2];
  }

  return 0.5F * (currents[n / 2 - 1] + currents[n / 2]);
}
