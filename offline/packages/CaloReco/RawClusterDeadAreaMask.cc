#include "RawClusterDeadAreaMask.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDeadMap.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>  // for pair
#include <vector>

RawClusterDeadAreaMask::RawClusterDeadAreaMask(const std::string &name)
  : SubsysReco(name)
  , m_detector("NONE")
  , m_deadTowerMaskHalfWidth(1.6)
  , m_rawClusters(nullptr)
  , m_deadMap(nullptr)
  , m_calibTowers(nullptr)
  , m_geometry(nullptr)
{
}

int RawClusterDeadAreaMask::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterDeadAreaMask::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << Name() << "::" << m_detector << "::process_event - Entry" << std::endl;
  }
  int nMasked = 0;

  const int eta_bins = m_geometry->get_etabins();
  const int phi_bins = m_geometry->get_phibins();
  assert(eta_bins > 0);
  assert(phi_bins > 0);

  //loop over the clusters
  RawClusterContainer::ConstRange begin_end = m_rawClusters->getClusters();
  RawClusterContainer::ConstIterator iter;

  for (iter = begin_end.first; iter != begin_end.second;)
  {
    //    RawClusterDefs::keytype key = iter->first;
    const RawCluster *cluster = iter->second;

    std::vector<float> toweretas;
    std::vector<float> towerphis;
    std::vector<float> towerenergies;

    //loop over the towers in the cluster
    RawCluster::TowerConstRange towers = cluster->get_towers();
    RawCluster::TowerConstIterator toweriter;

    for (toweriter = towers.first;
         toweriter != towers.second;
         ++toweriter)
    {
      RawTower *tower = m_calibTowers->getTower(toweriter->first);
      assert(tower);

      int towereta = tower->get_bineta();
      int towerphi = tower->get_binphi();
      double towerenergy = tower->get_energy();

      //put the etabin, phibin, and energy into the corresponding vectors
      toweretas.push_back(towereta);
      towerphis.push_back(towerphi);
      towerenergies.push_back(towerenergy);
    }

    int ntowers = toweretas.size();
    assert(ntowers >= 1);

    //loop over the towers to determine the energy
    //weighted eta and phi position of the cluster

    float etamult = 0;
    float etasum = 0;
    float phimult = 0;
    float phisum = 0;

    for (int j = 0; j < ntowers; j++)
    {
      float energymult = towerenergies.at(j) * toweretas.at(j);
      etamult += energymult;
      etasum += towerenergies.at(j);

      int phibin = towerphis.at(j);

      if (phibin - towerphis.at(0) < -phi_bins / 2.0)
        phibin += phi_bins;
      else if (phibin - towerphis.at(0) > +phi_bins / 2.0)
        phibin -= phi_bins;
      assert(std::abs(phibin - towerphis.at(0)) <= phi_bins / 2.0);

      energymult = towerenergies.at(j) * phibin;
      phimult += energymult;
      phisum += towerenergies.at(j);
    }

    float avgphi = phimult / phisum;
    float avgeta = etamult / etasum;

    if (Verbosity() > VERBOSITY_MORE)
    {
      std::cout << Name() << "::" << m_detector << "::process_event - "
                << "process cluster at average location " << avgeta << "," << avgphi << " : ";
      cluster->identify();
    }

    // rejection if close to a dead tower
    bool rejecCluster = false;

    for (int search_eta = ceil(avgeta - m_deadTowerMaskHalfWidth); search_eta <= floor(avgeta + m_deadTowerMaskHalfWidth); ++search_eta)
    {
      for (int search_phi = ceil(avgphi - m_deadTowerMaskHalfWidth); search_phi <= floor(avgphi + m_deadTowerMaskHalfWidth); ++search_phi)
      {
        int ieta = search_eta;
        int iphi = search_phi;

        if (ieta >= eta_bins) continue;
        if (ieta < 0) continue;

        if (iphi >= phi_bins) iphi -= phi_bins;
        if (iphi < 0) iphi += phi_bins;

        const bool isDead = m_deadMap->isDeadTower(ieta, iphi);

        // dead tower found in cluster
        if (Verbosity() > VERBOSITY_MORE)
        {
          std::cout << "\t"
                    << "tower " << ieta
                    << "," << iphi
                    << (isDead ? ": is dead." : "OK")
                    << std::endl;
        }

        if (isDead)
        {
          rejecCluster = true;
          break;
        }
      }

      if (rejecCluster) break;
    }

    //container operation
    if (rejecCluster)
    {
      if (Verbosity() > VERBOSITY_MORE)
      {
        std::cout << Name() << "::" << m_detector << "::process_event - "
                  << "reject cluster " << cluster->get_id() << std::endl;
        cluster->identify();
      }

      ++nMasked;
      m_rawClusters->getClustersMap().erase(iter++);
    }
    else
    {
      if (Verbosity() > VERBOSITY_MORE)
      {
        std::cout << Name() << "::" << m_detector << "::process_event - "
                  << "keep cluster " << cluster->get_id() << std::endl;
      }
      ++iter;
    }

  }  //  for (iter = begin_end.first; iter != begin_end.second;)

  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << Name() << "::" << m_detector << "::process_event - masked "
              << nMasked << " clusters. Final cluster containers has "
              << m_rawClusters->size() << " clusters"
              << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterDeadAreaMask::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  //Get the DST Node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  //Check that it is there
  if (!dstNode)
  {
    std::cout << "DST Node missing, quitting" << std::endl;
    throw std::runtime_error("failed to find DST node in RawClusterDeadAreaMask::CreateNodeTree");
  }

  m_rawClusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_" + m_detector);
  if (!m_rawClusters)
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "CreateNodeTree "
                 "No "
              << m_detector << " Cluster Container found while in RawClusterDeadAreaMask, can't proceed!!!" << std::endl;
    topNode->print();
    throw std::runtime_error("failed to find CLUSTER node in RawClusterDeadAreaMask::CreateNodeTree");
  }

  m_calibTowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_" + m_detector);
  if (!m_calibTowers)
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "CreateNodeTree "
              << "No calibrated " << m_detector << " tower info found while in RawClusterDeadAreaMask, can't proceed!!!" << std::endl;
    throw std::runtime_error("failed to find TOWER_CALIB node in RawClusterDeadAreaMask::CreateNodeTree");
  }

  std::string towergeomnodename = "TOWERGEOM_" + m_detector;
  m_geometry = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!m_geometry)
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "CreateNodeTree"
              << ": Could not find node " << towergeomnodename << std::endl;
    throw std::runtime_error("failed to find TOWERGEOM node in RawClusterDeadAreaMask::CreateNodeTree");
  }

  const std::string deadMapName = "DEADMAP_" + m_detector;
  m_deadMap = findNode::getClass<RawTowerDeadMap>(topNode, deadMapName);
  if (m_deadMap)
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "CreateNodeTree"
              << " use dead map: ";
    m_deadMap->identify();
  }
}

int RawClusterDeadAreaMask::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
