#include "RawTowerDeadTowerInterp.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDeadMap.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv1.h>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <cassert>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <set>  // for _Rb_tree_const_iterator
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

RawTowerDeadTowerInterp::RawTowerDeadTowerInterp(const std::string &name)
  : SubsysReco(name)
{
}

int RawTowerDeadTowerInterp::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode",
                                                           "DST"));
  if (!dstNode)
  {
    std::cout << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerDeadTowerInterp::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity())
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "process_event"
              << "Process event entered" << std::endl;
  }

  double recovery_energy = 0;
  int recoverTower = 0;
  int deadTowerCnt = 0;

  if (m_deadTowerMap)
  {
    const int eta_bins = m_geometry->get_etabins();
    const int phi_bins = m_geometry->get_phibins();

    // assume cylindrical calorimeters for now. Need to add swithc for planary calorimeter
    assert(eta_bins > 0);
    assert(phi_bins > 0);

    for (const RawTowerDefs::keytype &key : m_deadTowerMap->getDeadTowers())
    {
      ++deadTowerCnt;

      if (Verbosity() >= VERBOSITY_MORE)
      {
        std::cout << Name() << "::" << m_detector << "::"
                  << "process_event"
                  << " - processing tower " << key;
      }

      // check deadmap validity
      RawTowerGeom *towerGeom = m_geometry->get_tower_geometry(key);
      if (towerGeom == nullptr)
      {
        std::cout << Name() << "::" << m_detector << "::"
                  << "process_event"
                  << " - invalid dead tower ID " << key << std::endl;

        exit(2);
      }

      const int bineta = towerGeom->get_bineta();
      const int binphi = towerGeom->get_binphi();

      if (Verbosity() >= VERBOSITY_MORE)
      {
        std::cout << " bin " << bineta << "-" << binphi;
        std::cout << ". Add neighbors: ";
      }

      assert(bineta >= 0);
      assert(bineta <= eta_bins);
      assert(binphi >= 0);
      assert(binphi <= phi_bins);

      // eight neighbors
      static const std::vector<std::pair<int, int>> neighborIndexs =
          {{+1, 0}, {+1, +1}, {0, +1}, {-1, +1}, {-1, 0}, {-1, -1}, {0, -1}, {+1, -1}};

      int n_neighbor = 0;
      double E_SumNeighbor = 0;
      for (const std::pair<int, int> &neighborIndex : neighborIndexs)
      {
        int ieta = bineta + neighborIndex.first;
        int iphi = binphi + neighborIndex.second;

        if (ieta >= eta_bins) continue;
        if (ieta < 0) continue;

        if (iphi >= phi_bins) iphi -= phi_bins;
        if (iphi < 0) iphi += phi_bins;

        assert(ieta >= 0);
        assert(ieta <= eta_bins);
        assert(iphi >= 0);
        assert(iphi <= phi_bins);

        if (m_deadTowerMap->isDeadTower(ieta, iphi)) continue;

        ++n_neighbor;

        assert(m_calibTowers);
        RawTower *neighTower =
            m_calibTowers->getTower(ieta, iphi);
        if (neighTower == nullptr) continue;

        if (Verbosity() >= VERBOSITY_MORE)
        {
          std::cout << neighTower->get_energy() << " (" << ieta << "-" << iphi << "), ";
        }

        E_SumNeighbor += neighTower->get_energy();

      }  // for (const pair<int, int> &neighborIndex : neighborIndexs)

      if (n_neighbor > 0 and E_SumNeighbor != 0)
      {
        RawTower *deadTower = m_calibTowers->getTower(key);

        if (deadTower == nullptr)
        {
          deadTower = new RawTowerv1();
        }
        assert(deadTower);

        deadTower->set_energy(E_SumNeighbor / n_neighbor);
        m_calibTowers->AddTower(key, deadTower);

        recovery_energy += deadTower->get_energy();
        ++recoverTower;

        if (Verbosity() >= VERBOSITY_MORE)
        {
          std::cout << " -> " << deadTower->get_energy() << " GeV @ " << deadTower->get_id();
        }

      }  // if (n_neighbor>0)
      else
      {
        if (Verbosity() >= VERBOSITY_MORE)
        {
          std::cout << "No neighbor towers found.";
        }

      }  // if (n_neighbor>0) -else

      if (Verbosity() >= VERBOSITY_MORE)
      {
        std::cout << std::endl;
      }
    }

  }  //if (m_deadTowerMap)
  else
  {
    static bool once = true;

    if (Verbosity() and once)
    {
      once = false;

      std::cout << Name() << "::" << m_detector << "::"
                << "process_event"
                << " - missing dead map node. Do nothing ..."
                << std::endl;
    }
  }

  if (Verbosity())
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "process_event"
              << "recovery_energy = " << recovery_energy
              << " GeV from " << recoverTower << " towers out of total " << deadTowerCnt << " dead towers"
              << ", output sum energy = "
              << m_calibTowers->getTotalEdep() << " GeV" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerDeadTowerInterp::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawTowerDeadTowerInterp::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "CreateNodes"
              << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find Run node in RawTowerDeadTowerInterp::CreateNodes");
  }

  const std::string deadMapName = "DEADMAP_" + m_detector;
  m_deadTowerMap = findNode::getClass<RawTowerDeadMap>(topNode, deadMapName);
  if (m_deadTowerMap)
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "CreateNodes"
              << " use dead map: ";
    m_deadTowerMap->identify();
  }

  const std::string geometry_node = "TOWERGEOM_" + m_detector;
  m_geometry = findNode::getClass<RawTowerGeomContainer>(topNode, geometry_node);
  if (!m_geometry)
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "CreateNodes"
              << " " << geometry_node << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + geometry_node + " node in RawTowerDeadTowerInterp::CreateNodes");
  }

  if (Verbosity() >= 1)
  {
    m_geometry->identify();
  }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << Name() << "::" << m_detector << "::"
              << "CreateNodes"
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in RawTowerDeadTowerInterp::CreateNodes");
  }

  const std::string rawTowerNodeName = "TOWER_" + _calib_tower_node_prefix + "_" + m_detector;
  m_calibTowers = findNode::getClass<RawTowerContainer>(dstNode, rawTowerNodeName);
  if (!m_calibTowers)
  {
    std::cout << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << " " << rawTowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + rawTowerNodeName + " node in RawTowerDeadTowerInterp::CreateNodes");
  }

  return;
}
