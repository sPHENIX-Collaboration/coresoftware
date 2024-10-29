
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
  // here we are basically going to see if the max cluster pT or max jet pT is above a certain threshold, otherwise we will abort the event
  // jet loop
  JetContainer *jets = findNode::getClass<JetContainer>(topNode, m_JetNodeName);
  if (!jets)
  {
    std::cout << "JetDSTSkimmer::process_event - Error - Can't find Jet Node " << m_JetNodeName << " therefore no selection can be made" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  float maxjetpt = 0;
  for (auto jet : *jets)
  {
    if (jet->get_pt() > maxjetpt)
    {
      maxjetpt = jet->get_pt();
    }
  }

  RawClusterContainer *_clusters = findNode::getClass<RawClusterContainer>(topNode, m_ClusterNodeName);
  if (!_clusters)
  {
    std::cout << "JetDSTSkimmer::process_event - Error - Can't find Cluster Node " << m_ClusterNodeName << " therefore no selection can be made" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  float maxclusterpt = 0;
  RawClusterContainer::Map clusterMap = _clusters->getClustersMap();
  for (auto &clusterPair : clusterMap)
  {
    RawCluster *recoCluster = (clusterPair.second);
    if (recoCluster->get_energy() > maxclusterpt)
    {
      maxclusterpt = recoCluster->get_energy();
    }
  }
  // place holder for identifying background events
  bool isBackground = isBackgroundEvent();

  bool keepEvent = false;
  if (maxjetpt > m_minJetPt || maxclusterpt > m_minClusterPt)
  {
    keepEvent = true;
  }
  else
  {
    // verbose output, verbose 2 also will show the event abortion by this module
    if (Verbosity() > 1)
    {
      std::cout << "JetDSTSkimmer::process_event - Event Rejected - Max Jet pT: " << maxjetpt << " Max Cluster pT: " << maxclusterpt << std::endl;
    }
  }
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
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool JetDSTSkimmer::isBackgroundEvent()
{
  // place holder for identifying background events
  return false;
}
