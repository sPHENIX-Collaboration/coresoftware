#include "EpdReco.h"

#include "EpdGeomV2.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfoDefs.h>

#include <cdbobjects/CDBTTree.h> // for CDBTTree

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h> // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h> // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h> // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h> // for PHWHERE

#include <TSystem.h>

#include <array> // for array
#include <cstdlib> // for exit
#include <iostream>

EpdReco::EpdReco(const std::string &name)
  : SubsysReco(name)
{
  FillTilePhi0Array();
  FillTilePhiArray();
}

int EpdReco::InitRun(PHCompositeNode *topNode) {

  if (!m_overrideCalibName) {
    m_calibName = "SEPD_NMIP_CALIB";
  }
  if (!m_overrideFieldName) {
    m_fieldname = "sepd_calib";
  }
  std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
  if (!calibdir.empty()) {
    cdbttree = new CDBTTree(calibdir);
  } else {
    std::cout << "EpdReco::::InitRun No calibration file for domain "
              << m_calibName << " found" << std::endl;
    exit(1);
  }

  CreateNodes(topNode);

  PHNodeIterator node_itr(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(
      node_itr.findFirst("PHCompositeNode", "RUN"));
  if (!runNode) {
    std::cout << PHWHERE << "RUN Node not found - that is fatal" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  PHNodeIterator runiter(runNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(
      runiter.findFirst("PHCompositeNode", m_Detector));
  if (!DetNode) {
    DetNode = new PHCompositeNode(m_Detector);
    runNode->addNode(DetNode);
  }

  EpdGeom *epdGeom = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
  if (!epdGeom) {
    epdGeom = new EpdGeomV2();
    PHIODataNode<PHObject> *newNode =
        new PHIODataNode<PHObject>(epdGeom, "TOWERGEOM_EPD", "PHObject");
    DetNode->addNode(newNode);
  }

  // fill epd geometry
  unsigned int epdchannels = 744;
  for (unsigned int ch = 0; ch < epdchannels; ch++) {
    unsigned int thiskey = TowerInfoDefs::encode_epd(ch);
    epdGeom->set_z(thiskey, GetTileZ(TowerInfoDefs::get_epd_arm(thiskey)));
    epdGeom->set_r(thiskey, GetTileR(TowerInfoDefs::get_epd_rbin(thiskey)));
    if (TowerInfoDefs::get_epd_rbin(thiskey) == 0) {
      epdGeom->set_phi0(thiskey,
                        GetTilePhi0(TowerInfoDefs::get_epd_phibin(thiskey)));
    } else {
      epdGeom->set_phi(thiskey,
                       GetTilePhi(TowerInfoDefs::get_epd_phibin(thiskey)));
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int EpdReco::process_event(PHCompositeNode *topNode) {
  if (Verbosity() > 1) {
    std::cout << "EpdReco::process_event -- entered" << std::endl;
  }

  TowerInfoContainer *_sepd_towerinfo =
      findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_SEPD");
  unsigned int ntowers = 0;
  if (_sepd_towerinfo) {
    ntowers = _sepd_towerinfo->size();
  }
  if (ntowers != 744) {
    std::cout << "sEPD container has unexpected size - exiting now!"
              << std::endl;
    exit(1);
  }

  TowerInfoContainer *m_TowerInfoContainer_calib =
      findNode::getClass<TowerInfoContainer>(topNode,
                                             m_TowerInfoNodeName_calib);
  if (!m_TowerInfoContainer_calib) {
    std::cout << PHWHERE << "Could not locate TowerInfoContainer node "
              << m_TowerInfoNodeName_calib << std::endl;
    exit(1);
  }

  for (unsigned int ch = 0; ch < ntowers; ch++) {
    float ch_time = _sepd_towerinfo->get_tower_at_channel(ch)->get_time_float();
    float ch_adc = _sepd_towerinfo->get_tower_at_channel(ch)->get_energy();
    float ch_mpv = cdbttree->GetFloatValue(ch, m_fieldname);
    double ch_nmip = ch_adc / ch_mpv;
    m_TowerInfoContainer_calib->get_tower_at_channel(ch)->set_energy(ch_nmip);
    m_TowerInfoContainer_calib->get_tower_at_channel(ch)->set_time_float(
        ch_time);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

float EpdReco::GetTilePhi(int thisphi) {
  return tilephi.at(thisphi);
}

float EpdReco::GetTilePhi0(int thisphi0) {
  return tilephi0.at(thisphi0);
}

float EpdReco::GetTileR(int thisr) {
  static const std::array<float,16> tileR {
      6.8,    11.2,   15.6,   20.565, 26.095, 31.625, 37.155, 42.685,
      48.215, 53.745, 59.275, 64.805, 70.335, 75.865, 81.395, 86.925};
  return tileR.at(thisr);
}

float EpdReco::GetTileZ(int thisz) {
  static const std::array<float,2> tileZ {-316.0, 316.0};
  return tileZ.at(thisz);
}

void EpdReco::CreateNodes(PHCompositeNode *topNode) {
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode =
      dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(
      dstiter.findFirst("PHCompositeNode", m_Detector));
  if (!DetNode) {
    DetNode = new PHCompositeNode(m_Detector);
    dstNode->addNode(DetNode);
  }

  TowerInfoContainer *m_TowerInfoContainer_calib =
      findNode::getClass<TowerInfoContainer>(DetNode,
                                             m_TowerInfoNodeName_calib);
  if (m_TowerInfoContainer_calib == nullptr) {
    m_TowerInfoContainer_calib =
        new TowerInfoContainerv4(TowerInfoContainer::DETECTOR::SEPD);
    PHIODataNode<PHObject> *TowerInfoNodecalib = new PHIODataNode<PHObject>(
        m_TowerInfoContainer_calib, m_TowerInfoNodeName_calib, "PHObject");
    DetNode->addNode(TowerInfoNodecalib);
  }

  return;
}

void EpdReco::FillTilePhiArray()
{
  size_t i = 0;
  for (float & iter : tilephi)
  {
    iter = ((2 * i + 1) * M_PI) / 24;
    i++;
  }
  return;
}

void EpdReco::FillTilePhi0Array()
{
  size_t i = 0;
  for (float & iter : tilephi0)
  {
    iter = ((2 * i + 1) * M_PI) / 12;
    i++;
  }
  return;
}
