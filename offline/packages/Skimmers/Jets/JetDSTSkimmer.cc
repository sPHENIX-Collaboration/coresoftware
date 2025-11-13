
#include "JetDSTSkimmer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

#include <phool/getClass.h>

#include <iostream>                        // for operator<<, basic_ostream
#include <map>                             // for operator!=, _Rb_tree_iterator
#include <utility>                         // for pair

//____________________________________________________________________________..
JetDSTSkimmer::JetDSTSkimmer(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int JetDSTSkimmer::process_event(PHCompositeNode *topNode)
{
  // If no thresholds are configured, keep all events
  if (m_JetNodePts.empty() && m_ClusterNodePts.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // Check all configured jet nodes; pass if any jet exceeds its node threshold
  bool passedJet = false;
  for (const auto &jetNodeAndThreshold : m_JetNodePts)
  {
    const std::string &jetNodeName = jetNodeAndThreshold.first;
    const float jetThreshold = jetNodeAndThreshold.second;

    JetContainer *jets = findNode::getClass<JetContainer>(topNode, jetNodeName);
    if (!jets)
    {
      if (Verbosity() > 0)
      {
        std::cout << "JetDSTSkimmer::process_event - Warning - Can't find Jet Node Not Skimming on this: " << jetNodeName << std::endl;
      }
      continue;
    }
    for (auto *jet : *jets)
    {
      if (jet->get_pt() > jetThreshold)
      {
        passedJet = true;
        break;
      }
    }
    if (passedJet)
    {
      break;
    }
  }

  // Check all configured cluster nodes; pass if any cluster exceeds its node threshold
  bool passedCluster = false;
  for (const auto &clusterNodeAndThreshold : m_ClusterNodePts)
  {
    const std::string &clusterNodeName = clusterNodeAndThreshold.first;
    const float clusterThreshold = clusterNodeAndThreshold.second;

    RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, clusterNodeName);
    if (!clusters)
    {
      if (Verbosity() > 0)
      {
        std::cout << "JetDSTSkimmer::process_event - Warning - Can't find Cluster Node " << clusterNodeName << std::endl;
      }
      continue;
    }
    RawClusterContainer::Map clusterMap = clusters->getClustersMap();
    for (auto &clusterPair : clusterMap)
    {
      RawCluster *recoCluster = (clusterPair.second);
      if (recoCluster->get_energy() > clusterThreshold)
      {
        passedCluster = true;
        break;
      }
    }
    if (passedCluster)
    {
      break;
    }
  }

  // place holder for identifying background events
  bool isBackground = isBackgroundEvent();

  bool keepEvent = (passedJet || passedCluster);
  if (isBackground)
  {
    keepEvent = false;
    if (Verbosity() > 1)
    {
      std::cout << "JetDSTSkimmer::process_event - Event Rejected - Background Event" << std::endl;
    }
  }

  if (!keepEvent)
  {
    if (Verbosity() > 1)
    {
      std::cout << "JetDSTSkimmer::process_event - Event Rejected - No jets or clusters above thresholds" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool JetDSTSkimmer::isBackgroundEvent()
{
  // place holder for identifying background events
  return false;
}
