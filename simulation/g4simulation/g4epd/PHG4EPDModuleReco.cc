
#include "PHG4EPDModuleReco.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>

#include <epd/EPDDefs.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for hit_idbits

#include <phparameter/PHParameterInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <map>      // for _Rb_tree_const_iterator
#include <utility>  // for pair

PHG4EPDModuleReco::PHG4EPDModuleReco(const std::string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();
}

int PHG4EPDModuleReco::InitRun(PHCompositeNode *topNode)
{
  tmin = get_double_param("tmin");
  tmax = get_double_param("tmax");
  m_DeltaT = get_double_param("delta_t");
  m_EpdMpv = get_double_param("epdmpv");

  CreateNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4EPDModuleReco::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, m_Hitnodename);
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << m_Hitnodename << std::endl;
    exit(1);
  }

  TowerInfoContainer *m_TowerInfoContainer = findNode::getClass<TowerInfoContainer>(topNode, m_TowerInfoNodeName);
  if (!m_TowerInfoContainer)
  {
    std::cout << PHWHERE << "Could not locate TowerInfoContainer node " << m_TowerInfoNodeName << std::endl;
    exit(1);
  }

  TowerInfoContainer *m_TowerInfoContainer_calib = findNode::getClass<TowerInfoContainer>(topNode, m_TowerInfoNodeName_calib);
  if (!m_TowerInfoContainer_calib)
  {
    std::cout << PHWHERE << "Could not locate TowerInfoContainer node " << m_TowerInfoNodeName_calib << std::endl;
    exit(1);
  }

  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    if (hiter->second->get_t(0) > tmax)
    {
      continue;
    }
    if (hiter->second->get_t(1) < tmin)
    {
      continue;
    }
    if (hiter->second->get_t(1) - hiter->second->get_t(0) > m_DeltaT)
    {
      continue;
    }

    int sim_tileid = (hiter->second->get_hit_id() >> PHG4HitDefs::hit_idbits);
    int this_tile = EPDDefs::get_tileid(sim_tileid);
    int this_sector = EPDDefs::get_sector(sim_tileid);
    int this_arm = EPDDefs::get_arm(sim_tileid);

    m_EpdTile_e[this_arm][this_sector][this_tile] += hiter->second->get_light_yield();
    m_EpdTile_Calib_e[this_arm][this_sector][this_tile] += hiter->second->get_light_yield() / m_EpdMpv;

  }  // end loop over g4hits

  for (unsigned int k = 0; k < 2; k++)
  {
    for (int i = 0; i < 12; i++)
    {
      for (int j = 0; j < 31; j++)
      {
        unsigned int globalphi = Getphimap(j) + 2 * i;
        unsigned int r = Getrmap(j);

        unsigned int key = globalphi + (r << 10U) + (k << 20U);

        unsigned int ch = m_TowerInfoContainer->decode_key(key);

        m_TowerInfoContainer->get_tower_at_channel(ch)->set_energy(m_EpdTile_e[k][i][j]);

        m_TowerInfoContainer_calib->get_tower_at_channel(ch)->set_energy(m_EpdTile_Calib_e[k][i][j]);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4EPDModuleReco::SetDefaultParameters()
{
  set_default_double_param("tmax", 60.0);
  set_default_double_param("tmin", -20.0);
  set_default_double_param("delta_t", 100.);
  set_default_double_param("epdmpv", 1.);
  return;
}

void PHG4EPDModuleReco::set_timing_window(const double tmi, const double tma)
{
  set_double_param("tmin", tmi);
  set_double_param("tmax", tma);
}

int PHG4EPDModuleReco::Getrmap(int rindex)
{
  int rmap[31] = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15};

  return rmap[rindex];
}

int PHG4EPDModuleReco::Getphimap(int phiindex)
{
  int phimap[31] = {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};

  return phimap[phiindex];
}

void PHG4EPDModuleReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_Detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_Detector);
    dstNode->addNode(DetNode);
  }

  m_TowerInfoNodeName = "TOWERINFO_" + m_EPDSimTowerNodePrefix + "_" + m_Detector; // detector name and prefix are set by now
  TowerInfoContainer *m_TowerInfoContainer = findNode::getClass<TowerInfoContainer>(DetNode, m_TowerInfoNodeName);
  if (m_TowerInfoContainer == nullptr)
  {
    m_TowerInfoContainer = new TowerInfoContainerv1(TowerInfoContainer::DETECTOR::SEPD);
    PHIODataNode<PHObject> *TowerInfoNode = new PHIODataNode<PHObject>(m_TowerInfoContainer, m_TowerInfoNodeName, "PHObject");
    DetNode->addNode(TowerInfoNode);
  }

  m_TowerInfoNodeName_calib = "TOWERINFO_" + m_EPDCalibTowerNodePrefix + "_" + m_Detector; // detector name and prefix are set by now
  TowerInfoContainer *m_TowerInfoContainer_calib = findNode::getClass<TowerInfoContainer>(DetNode, m_TowerInfoNodeName_calib);
  if (m_TowerInfoContainer_calib == nullptr)
  {
    m_TowerInfoContainer_calib = new TowerInfoContainerv1(TowerInfoContainer::DETECTOR::SEPD);
    PHIODataNode<PHObject> *TowerInfoNodecalib = new PHIODataNode<PHObject>(m_TowerInfoContainer_calib, m_TowerInfoNodeName_calib, "PHObject");
    DetNode->addNode(TowerInfoNodecalib);
  }

  return;
}
int PHG4EPDModuleReco::ResetEvent(PHCompositeNode * /*topNode*/)
{
  // this only works for initializing to zero
  m_EpdTile_Calib_e = {};
  m_EpdTile_e = {};
  // if you ever want to initialize to a different value, do it this way:
  //   for (auto &entry : epd_tile_calib_e)
  //   {
  //     for (auto &entry1 : entry)
  //     {
  // entry1.fill(NAN);
  //     }
  //   }
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4EPDModuleReco::Detector(const std::string &detector)
{
  m_Detector = detector;
  m_Hitnodename = "G4HIT_" + m_Detector;
}
