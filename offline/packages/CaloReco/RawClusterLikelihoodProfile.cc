#include "RawClusterLikelihoodProfile.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>

#include <iostream>
#include <vector>

RawClusterLikelihoodProfile::RawClusterLikelihoodProfile(const std::string &name)
  : SubsysReco(name)
{
}

RawClusterLikelihoodProfile::~RawClusterLikelihoodProfile() = default;

int RawClusterLikelihoodProfile::Init(PHCompositeNode *topNode)
{
  // load profile for prob calculation
  cdfcalc = new ClusterCDFCalculator();
  cdfcalc->LoadProfile(m_profile_name);

  if (m_inputNodeName == m_outputNodeName)
  {
    std::cout << "RawClusterLikelihoodProfile::Init: Same inputNodeName and outputNodeName, putting inplace to true for overwriting" << std::endl;
    inplace = true;
  }
  CreateNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterLikelihoodProfile::process_event(PHCompositeNode *topNode)
{
  std::string clusterNodeName = m_inputNodeName;
  RawClusterContainer *clusterContainer = findNode::getClass<RawClusterContainer>(topNode, clusterNodeName);
  if (!clusterContainer)
  {
    std::cout << "RawClusterLikelihoodProfile::process_event::Could not locate input cluster node " << clusterNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (inplace)
  {
    _clusters = clusterContainer;
  }
  else
  {
    _clusters->Reset();
    RawClusterContainer::Map clusterMap = clusterContainer->getClustersMap();
    for (auto &clusterPair : clusterMap)
    {
      RawCluster *cluster = (clusterPair.second);
      float clusterE = cluster->get_energy();
      if (clusterE < m_min_cluster_e)
      {
        continue;
      }
      RawCluster *profileCluster = (RawCluster *) cluster->CloneMe();
      _clusters->AddCluster(profileCluster);
    }
  }

  std::string towerNodeName = m_towerNodeName;
  TowerInfoContainer *emcTowerContainer = findNode::getClass<TowerInfoContainer>(topNode, towerNodeName);
  if (!emcTowerContainer)
  {
    std::cout << "RawClusterLikelihoodProfile::process_event Could not locate tower node " << towerNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  RawClusterContainer::Map clusterMap = _clusters->getClustersMap();
  for (auto &clusterPair : clusterMap)
  {
    RawCluster *cluster = clusterPair.second;
    cluster->set_prob(-1);
    CLHEP::Hep3Vector vertex(0, 0, 0);
    CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*cluster, vertex);
    float ET = E_vec_cluster_Full.perp();
    if (ET < minET)
    {
      continue;
    }
    const RawCluster::TowerMap tower_map = cluster->get_towermap();

    std::vector<float> shower_shapes = cluster->get_shower_shapes(m_tower_thres_e);
    int ieta_center_of_gravity = std::floor(shower_shapes[4] + 0.5);
    int iphi_center_of_gravity = std::floor(shower_shapes[5] + 0.5);

    std::vector<double> input;
    int vectorSize = inputDimx * inputDimy;
    input.resize(vectorSize, 0);

    int xlength = int((inputDimx - 1) / 2);
    int ylength = int((inputDimy - 1) / 2);
    if (ieta_center_of_gravity - ylength < 0 || ieta_center_of_gravity + ylength >= 96)
    {
      continue;
    }
    for (int ieta = ieta_center_of_gravity - ylength; ieta <= ieta_center_of_gravity + ylength; ieta++)
    {
      for (int iphi = iphi_center_of_gravity - xlength; iphi <= iphi_center_of_gravity + xlength; iphi++)
      {
        int mappediphi = iphi;

        if (mappediphi < 0)
        {
          mappediphi += 256;
        }
        if (mappediphi > 255)
        {
          mappediphi -= 256;
        }
        unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, mappediphi);
        TowerInfo *towerinfo = emcTowerContainer->get_tower_at_key(towerinfokey);
        if (!towerinfo)
        {
          std::cout << "No towerinfo for tower key " << towerinfokey << std::endl;
          std::cout << "ieta: " << ieta << " iphi: " << mappediphi << std::endl;
          continue;
        }
        int index = (ieta - ieta_center_of_gravity + ylength) * inputDimx + iphi - iphi_center_of_gravity + xlength;
        input.at(index) = towerinfo->get_energy();
      }
    }
    double prob = cdfcalc->GetCDF(input, m_profile_dimension);
    cluster->set_prob(prob);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterLikelihoodProfile::End(PHCompositeNode * /*topNode*/)
{
  delete cdfcalc;
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterLikelihoodProfile::Print(const std::string &what) const
{
  std::cout << "RawClusterLikelihoodProfile::Print(const std::string &what) const Printing info for " << what << std::endl;
  return;
}

void RawClusterLikelihoodProfile::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Grab the CEMC node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
  }

  // Get the _det_name subnode
  PHCompositeNode *cemcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "CEMC"));

  // Check that it is there
  if (!cemcNode)
  {
    cemcNode = new PHCompositeNode("CEMC");
    dstNode->addNode(cemcNode);
  }
  std::string clusterNodeName = m_outputNodeName;
  _clusters = findNode::getClass<RawClusterContainer>(dstNode, clusterNodeName);
  if (!_clusters)
  {
    _clusters = new RawClusterContainer();
    PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, clusterNodeName, "PHObject");
    cemcNode->addNode(clusterNode);
  }
}
