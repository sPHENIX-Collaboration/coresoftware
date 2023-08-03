#include "TowerInfoDeadHotMask.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerDeadMap.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <fstream>
#include <iostream>
#include <set>                               // for _Rb_tree_const_iterator
#include <sstream>
#include <stdexcept>
#include <string>

TowerInfoDeadHotMask::TowerInfoDeadHotMask(const std::string &name)
  : SubsysReco(name)
  , m_detector("NONE")
  , m_deadMap(nullptr)
  , m_calibTowerInfos(nullptr)
  , m_geometry(nullptr)
{
}

int TowerInfoDeadHotMask::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int TowerInfoDeadHotMask::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() >= VERBOSITY_SOME)
  {
    std::cout << Name() << "::" << m_detector << "::process_event - Entry" << std::endl;
  }

  RawTowerDeadMap::Map map = m_deadMap->getDeadTowers();

  RawTowerDeadMap::Map::iterator itr;

  for(itr = map.begin(); itr != map.end(); ++itr)
  {
    int iphi = m_geometry->get_tower_geometry(*itr)->get_binphi();
    int ieta = m_geometry->get_tower_geometry(*itr)->get_bineta();
    unsigned int key = TowerInfoDefs::encode_emcal(ieta, iphi);
    m_calibTowerInfos->get_tower_at_key(key)->set_energy(0.);
    m_calibTowerInfos->get_tower_at_key(key)->set_time(-10);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TowerInfoDeadHotMask::CreateNodeTree(PHCompositeNode *topNode)
{
  m_calibTowerInfos = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_" + m_detector);
  if (!m_calibTowerInfos)
  {
    std::cout << Name() << "::" << m_detector << "::" << "CreateNodeTree " << "No calibrated " << m_detector << " tower info found while in TowerInfoDeadHotMask, can't proceed!!!" << std::endl;
    throw std::runtime_error("failed to find TOWERINFO_CALIB node in TowerInfoDeadHotMask::CreateNodeTree");
  }

  std::string towergeomnodename = "TOWERGEOM_" + m_detector;
  m_geometry = findNode::getClass<RawTowerGeomContainer>(topNode, towergeomnodename);
  if (!m_geometry)
  {
    std::cout << Name() << "::" << m_detector << "::" << "CreateNodeTree" << ": Could not find node " << towergeomnodename << std::endl;
    throw std::runtime_error("failed to find TOWERGEOM node in TowerInfoDeadHotMask::CreateNodeTree");
  }

  const std::string deadMapName = "DEADMAP_" + m_detector;
  m_deadMap = findNode::getClass<RawTowerDeadMap>(topNode, deadMapName);
  if (m_deadMap)
  {
    std::cout << Name() << "::" << m_detector << "::" << "CreateNodeTree" << " use dead map: ";
    m_deadMap->identify();
  }
}

int TowerInfoDeadHotMask::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
