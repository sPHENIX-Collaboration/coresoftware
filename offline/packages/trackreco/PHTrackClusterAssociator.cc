
#include "PHTrackClusterAssociator.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackCaloClusterMap_v1.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <phgeom/PHGeomUtility.h>

namespace
{
  template <class T>
  inline constexpr T deltaPhi(const T& dphi)
  {
    if (dphi > M_PI)
      return dphi - 2. * M_PI;
    else if (dphi <= -M_PI)
      return dphi + 2. * M_PI;
    else
      return dphi;
  }
}  // namespace

//____________________________________________________________________________..
PHTrackClusterAssociator::PHTrackClusterAssociator(const std::string& name)
  : SubsysReco(name)
{
  m_caloNames.push_back("CEMC");
  m_caloNames.push_back("HCALIN");
  m_caloNames.push_back("HCALOUT");
}

//____________________________________________________________________________..
PHTrackClusterAssociator::~PHTrackClusterAssociator()
{
}

//____________________________________________________________________________..
int PHTrackClusterAssociator::InitRun(PHCompositeNode* topNode)
{
  int ret = getNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  ret = createNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  return ret;
}

//____________________________________________________________________________..
int PHTrackClusterAssociator::process_event(PHCompositeNode* topNode)
{
  if (!m_calosAvailable)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  for (int layer = 0; layer < m_nCaloLayers; layer++)
  {
    if (Verbosity() > 1)
    {
      std::cout << "Processing calo layer "
                << m_caloNames.at(layer) << std::endl;
    }

    int ret = matchTracks(topNode, layer);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  if (Verbosity() > 3)
  {
    for (const auto& [track, clustervec] : *m_trackClusterMap)
    {
      track->identify();
      std::cout << " has clusters associated to it : " << std::endl;
      for (auto cluster : clustervec)
      {
        cluster->identify();
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackClusterAssociator::matchTracks(PHCompositeNode* topNode,
                                          const int caloLayer)
{
  if (getCaloNodes(topNode, caloLayer) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /// Default to using calo radius
  double caloRadius = m_towerGeomContainer->get_radius();
  if (m_caloRadii.find(m_caloNames.at(caloLayer)) != m_caloRadii.end())
  {
    caloRadius = m_caloRadii.find(m_caloNames.at(caloLayer))->second;
  }

  for (const auto& [key, track] : *m_trackMap)
  {
    const SvtxTrackState* state = track->get_state(caloRadius);
    const float statex = state->get_x();
    const float statey = state->get_y();
    const float statez = state->get_z();
    const float statephi = atan2(statey, statex);
    const float stateeta = asinh(statez /
                                 sqrt(statex * statex + statey * statey));

    const int vertexid = track->get_vertex_id();
    const auto vertex = m_vertexMap->get(vertexid);
    const auto cluster = getCluster(statephi, stateeta, vertex);
    if (Verbosity() > 1)
    {
      if (!cluster)
      {
        std::cout << "no cluster found, continuing to next track"
                  << std::endl;
        continue;
      }
      else
      {
        std::cout << "matching cluster " << cluster->get_id() << " to track " << track->get_id() << std::endl;
      }
    }

    m_trackClusterMap->insert(track, cluster);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

RawCluster* PHTrackClusterAssociator::getCluster(double phi,
                                                 double eta,
                                                 SvtxVertex* vertex)
{
  double minR = std::numeric_limits<double>::max();
  auto clusterMap = m_clusterContainer->getClustersMap();
  RawCluster* returncluster = nullptr;
  Acts::Vector3 vert = Acts::Vector3::Zero();
  if (vertex)
  {
    vert(0) = vertex->get_x();
    vert(1) = vertex->get_y();
    vert(2) = vertex->get_z();
  }

  for (const auto& [key, cluster] : clusterMap)
  {
    const auto clusterEta =
        RawClusterUtility::GetPseudorapidity(*cluster,
                                             CLHEP::Hep3Vector(vert(0), vert(1), vert(2)));
    const auto dphi = deltaPhi(phi - cluster->get_phi());
    const auto deta = eta - clusterEta;
    const auto r = sqrt(pow(dphi, 2) + pow(deta, 2));

    if (r < minR)
    {
      minR = r;
      returncluster = cluster;
    }
  }

  return returncluster;
}

int PHTrackClusterAssociator::getNodes(PHCompositeNode* topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_trackMap)
  {
    std::cout << PHWHERE << "No track map, can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
  {
    std::cout << PHWHERE << "No vertex map, can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackClusterAssociator::getCaloNodes(PHCompositeNode* topNode,
                                           const int caloLayer)
{
  std::string towerGeoNodeName = "TOWERGEOM_" + m_caloNames.at(caloLayer);
  std::string towerNodeName = "TOWER_CALIB_" + m_caloNames.at(caloLayer);
  std::string clusterNodeName = "CLUSTER_" + m_caloNames.at(caloLayer);

  m_towerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, towerGeoNodeName.c_str());

  m_towerContainer = findNode::getClass<RawTowerContainer>(topNode, towerNodeName.c_str());

  m_clusterContainer = findNode::getClass<RawClusterContainer>(topNode, clusterNodeName.c_str());

  if (m_useCemcPosRecalib and
      m_caloNames.at(caloLayer).compare("CEMC") == 0)
  {
    std::string nodeName = "CLUSTER_POS_COR_" + m_caloNames.at(caloLayer);
    m_clusterContainer = findNode::getClass<RawClusterContainer>(topNode, nodeName.c_str());
  }

  if (!m_towerGeomContainer or !m_towerContainer or !m_clusterContainer)
  {
    std::cout << PHWHERE
              << "Calo geometry and/or cluster container not found on node tree. Track-Calo cluster map won't be filled."
              << std::endl;
    m_calosAvailable = false;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackClusterAssociator::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in PHActsSiliconSeeding::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_trackClusterMap = findNode::getClass<SvtxTrackCaloClusterMap>(topNode, "TrackCaloClusterMap");
  if (!m_trackClusterMap)
  {
    m_trackClusterMap = new SvtxTrackCaloClusterMap_v1;
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(
        m_trackClusterMap, "TrackCaloClusterMap", "PHObject");
    svtxNode->addNode(node);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTrackClusterAssociator::ResetEvent(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int PHTrackClusterAssociator::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
