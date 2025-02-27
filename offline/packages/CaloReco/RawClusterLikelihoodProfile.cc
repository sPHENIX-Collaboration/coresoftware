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

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <ffamodules/CDBInterface.h>
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

int RawClusterLikelihoodProfile::Init(PHCompositeNode *topNode)
{
  cdfcalc = new ClusterCDFCalculator();
  cdfcalc->LoadProfile(m_profile_path);

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

  std::string globalvtxNodeName = m_globalvtxNodeName;
  GlobalVertexMap *globalvtxmap = findNode::getClass<GlobalVertexMap>(topNode,globalvtxNodeName);
  bool isglbvtx=true;
  if(!globalvtxmap || globalvtxmap->empty())
  {
    isglbvtx=false;
  }

  if(isglbvtx){
    GlobalVertex *bvertex= nullptr;
    for (GlobalVertexMap::ConstIter globaliter= globalvtxmap->begin(); globaliter != globalvtxmap->end(); ++globaliter)
    {
      bvertex = globaliter->second;
    }
    if(!bvertex){std::cout << "WARNING! RawClusterLikelihoodProfile::could not find globalvtxmap iter :: set vtx to (0,0,0)" << std::endl;}
    else if(bvertex){
      vz = bvertex->get_z();
      vy = bvertex->get_y();
      vx = bvertex->get_x();
    }
  }

  RawClusterContainer::Map clusterMap = _clusters->getClustersMap();
  for (auto &clusterPair : clusterMap)
  {
    RawCluster *cluster = clusterPair.second;
    cluster->set_prob(-1);
    cluster->set_merged_cluster_prob(-1);
    CLHEP::Hep3Vector vertex(vx, vy, vz);
    CLHEP::Hep3Vector E_vec_cluster_Full = RawClusterUtility::GetEVec(*cluster, vertex);
    double clusterE = E_vec_cluster_Full.mag();
    if (clusterE < m_min_cluster_e)
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
    std::pair<double, double> probpair = cdfcalc->GetCDF(input, clusterE, m_profile_dimension);
    double prob = probpair.first; 
    double prob_merged_cluster = probpair.second; 
    cluster->set_prob(prob);
    cluster->set_merged_cluster_prob(prob_merged_cluster);
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

  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in EmcRawTowerBuilder::CreateNodes");
  }

  PHCompositeNode *cemcNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "CEMC"));

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
