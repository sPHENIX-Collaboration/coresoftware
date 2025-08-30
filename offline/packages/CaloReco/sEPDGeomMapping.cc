#include "sEPDGeomMapping.h"

#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainerv1.h>
#include <calobase/RawTowerGeomv1.h>
#include <calobase/TowerInfoDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>
#include <iostream>

sEPDGeomMapping::sEPDGeomMapping(const std::string &name)
  : SubsysReco(name)
{
}

int sEPDGeomMapping::InitRun(PHCompositeNode *topNode)
{
  // Fill the phi arrays
  FillTilePhiArray();
  FillTilePhi0Array();

  CreateGeometry(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void sEPDGeomMapping::CreateGeometry(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << " Run Node missing, doing nothing." << std::endl;
    return;
  }

  PHNodeIterator runiter(runNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(runiter.findFirst("PHCompositeNode", "SEPD"));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode("SEPD");
    runNode->addNode(DetNode);
  }

  m_RawTowerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_SEPD");
  if (!m_RawTowerGeomContainer)
  {
    m_RawTowerGeomContainer = new RawTowerGeomContainerv1(RawTowerDefs::SEPD);
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(m_RawTowerGeomContainer, "TOWERGEOM_SEPD", "PHObject");
    DetNode->addNode(newNode);
  }
  else
  {
    std::cout << "sEPDGeomMapping::CreateGeometry - TOWERGEOM_SEPD already exists, not recreating" << std::endl;
    return;
  }


  unsigned int epdchannels = 744;
  for (unsigned int ch = 0; ch < epdchannels; ch++)
  {
    unsigned int towerkey = TowerInfoDefs::encode_epd(ch);
    int arm = TowerInfoDefs::get_epd_arm(towerkey);
    int rbin = TowerInfoDefs::get_epd_rbin(towerkey);
    int phibin = TowerInfoDefs::get_epd_phibin(towerkey);

    unsigned int key = TowerInfoDefs::get_sepd_geokey_at_channel(ch);

    float r = GetTileR(rbin);
    float phi = (rbin == 0) ? GetTilePhi0(phibin) : GetTilePhi(phibin);
    float z = GetTileZ(arm);

    float x = r * cos(phi);
    float y = r * sin(phi);

    RawTowerGeom *tg = new RawTowerGeomv1(key);
    tg->set_center_x(x);
    tg->set_center_y(y);
    tg->set_center_z(z);

    m_RawTowerGeomContainer->add_tower_geometry(tg);
  }

  if (Verbosity() > 0)
  {
    std::cout << "sEPDGeomMapping::CreateGeometry - created " << epdchannels 
              << " tower geometries for SEPD" << std::endl;
  }
}

void sEPDGeomMapping::FillTilePhiArray()
{
  tilephi = {{0.13089969, 0.39269908, 0.65449847, 0.91629786, 1.17809725, 1.43989663, 
              1.70169602, 1.96349541, 2.2252948, 2.48709418, 2.74889357, 3.01069296, 
              3.27249235, 3.53429174, 3.79609112, 4.05789051, 4.3196899, 4.58148929, 
              4.84328867, 5.10508806, 5.36688745, 5.62868684, 5.89048623, 6.15228561}};
}

void sEPDGeomMapping::FillTilePhi0Array()
{
  tilephi0 = {{0.26179939, 0.78539816, 1.30899694, 1.83259571, 2.35619449, 2.87979327, 
               3.40339204, 3.92699082, 4.45058959, 4.97418837, 5.49778714, 6.02138592}};
}

float sEPDGeomMapping::GetTilePhi(int thisphi) const
{
  if (thisphi < 0 || thisphi >= 24)
  {
    std::cout << "sEPDGeomMapping::GetTilePhi - invalid phibin " << thisphi << std::endl;
    return 0;
  }
  return tilephi[thisphi];
}

float sEPDGeomMapping::GetTilePhi0(int thisphi0) const
{
  if (thisphi0 < 0 || thisphi0 >= 12)
  {
    std::cout << "sEPDGeomMapping::GetTilePhi0 - invalid phibin " << thisphi0 << std::endl;
    return 0;
  }
  return tilephi0[thisphi0];
}

float sEPDGeomMapping::GetTileR(int thisr) const
{
  static const float tileR[16] = {6.8, 11.2, 15.6, 20.565, 26.095, 31.625, 37.155, 42.685, 
                                   48.215, 53.745, 59.275, 64.805, 70.335, 75.865, 81.395, 86.925};
  if (thisr < 0 || thisr >= 16)
  {
    std::cout << "sEPDGeomMapping::GetTileR - invalid rbin " << thisr << std::endl;
    return 0;
  }
  return tileR[thisr];
}

float sEPDGeomMapping::GetTileZ(int thisz) const
{
  static const float tileZ[2] = {-316.0, 316.0};
  if (thisz < 0 || thisz >= 2)
  {
    std::cout << "sEPDGeomMapping::GetTileZ - invalid arm " << thisz << std::endl;
    return 0;
  }
  return tileZ[thisz];
}

int sEPDGeomMapping::GetRmap(int tileindex) const
{
  static const int rmap[31] = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 
                                8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15};
  if (tileindex < 0 || tileindex >= 31)
  {
    std::cout << "sEPDGeomMapping::GetRmap - invalid tile index " << tileindex << std::endl;
    return 0;
  }
  return rmap[tileindex];
}

int sEPDGeomMapping::GetPhimap(int phiindex) const
{
  static const int phimap[31] = {0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 
                                  1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
  if (phiindex < 0 || phiindex >= 31)
  {
    std::cout << "sEPDGeomMapping::GetPhimap - invalid tile index " << phiindex << std::endl;
    return 0;
  }
  return phimap[phiindex];
}