#include "EventPlaneReco.h"

#include "Eventplaneinfo.h"
#include "EventplaneinfoMap.h"
#include "EventplaneinfoMapv1.h"
#include "Eventplaneinfov1.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <epd/EpdGeom.h>

#include <mbd/MbdGeom.h>
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <cdbobjects/CDBHistos.h>
#include <cdbobjects/CDBTTree.h>
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
#include <phool/recoConsts.h>

#include <array> // for array
#include <cfloat>
#include <cmath>
#include <cstdlib> // for exit
#include <iostream>
#include <set>     // for _Rb_tree_const_iterator
#include <utility> // for pair
#include <vector>  // for vector

EventPlaneReco::EventPlaneReco(const std::string &name) : SubsysReco(name) {
  south_q.resize(m_MaxOrder);
  north_q.resize(m_MaxOrder);

  for (auto &vec : south_q) {
    vec.resize(2);
  }

  for (auto &vec : north_q) {
    vec.resize(2);
  }
}

int EventPlaneReco::InitRun(PHCompositeNode *topNode) {
  if (!m_overrideSEPDMapName) {
    m_sEPDMapName = "SEPD_CHANNELMAP";
  }
  if (!m_overrideSEPDFieldName) {
    m_sEPDfieldname = "epd_channel_map";
  }
  std::string calibdir = CDBInterface::instance()->getUrl(m_sEPDMapName);
  if (!calibdir.empty()) {
    cdbttree = new CDBTTree(calibdir);
  } else {
    std::cout << "EventPlaneReco::::InitRun No sEPD mapping file for domain "
              << m_sEPDMapName << " found" << std::endl;
    exit(1);
  }
  vkey.clear();
  for (int i = 0; i < 768; i++) {

    int keymap = cdbttree->GetIntValue(i, m_sEPDfieldname);
    if (keymap == 999) {
      continue;
    }

    key = TowerInfoDefs::encode_epd(keymap);
    vkey.push_back(key);
  }

  return CreateNodes(topNode);
}

int EventPlaneReco::process_event(PHCompositeNode *topNode) {
  if (Verbosity() > 1) {
    std::cout << "EventPlaneReco::process_event -- entered" << std::endl;
  }

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  EventplaneinfoMap *epmap =
      findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!epmap) {
    std::cout << PHWHERE << "::ERROR - cannot find EventplaneinfoMap"
              << std::endl;
    exit(-1);
  }

  if (_sepdEpReco) {
    ResetMe();
    TowerInfoContainer *epd_towerinfo =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_SEPD");
    if (!epd_towerinfo) {
      std::cout << PHWHERE << "::ERROR - cannot find TOWERINFO_CALIB_SEPD"
                << std::endl;
      exit(-1);
    }
    EpdGeom *_epdgeom = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
    if (!_epdgeom) {
      std::cout << PHWHERE << "::ERROR - cannot find TOWERGEOM_EPD"
                << std::endl;
      exit(-1);
    }

    if (epd_towerinfo) {
      if (Verbosity()) {
        std::cout << "EventPlaneReco::process_event -  epd_towerinfo"
                  << std::endl;
      }

      unsigned int ntowers = epd_towerinfo->size();
      for (unsigned int ch = 0; ch < ntowers; ch++) {
        TowerInfo *_tower = epd_towerinfo->get_tower_at_channel(ch);
        float epd_e = _tower->get_energy();
        float epd_time = _tower->get_time_float();
        if (epd_time > 0.) // exclude ZS
        {
          if (epd_e < 0.5) // expecting Nmips
          {
            continue;
          }
          float tile_phi = _epdgeom->get_phi(vkey[ch]);
          int arm = TowerInfoDefs::get_epd_arm(vkey[ch]);
          float truncated_e =
              (epd_e < _epd_e) ? epd_e : _epd_e; // set cutoff at _epd_e
          if (arm == 0) {
            for (unsigned int order = 0; order < m_MaxOrder; order++) {
              double Cosine = cos(tile_phi * (double)(order + 1));
              double Sine = sin(tile_phi * (double)(order + 1));
              south_q[order][0] += truncated_e * Cosine; // south Qn,x
              south_q[order][1] += truncated_e * Sine;   // south Qn,y
            }
          } else if (arm == 1) {
            for (unsigned int order = 0; order < m_MaxOrder; order++) {
              double Cosine = cos(tile_phi * (double)(order + 1));
              double Sine = sin(tile_phi * (double)(order + 1));
              north_q[order][0] += truncated_e * Cosine; // north Qn,x
              north_q[order][1] += truncated_e * Sine;   // north Qn,y
            }
          }
        }
      }
    }
    for (unsigned int order = 0; order < m_MaxOrder; order++) {
      south_Qvec.emplace_back(south_q[order][0], south_q[order][1]);
      north_Qvec.emplace_back(north_q[order][0], north_q[order][1]);
    }

    if (epd_towerinfo) {
      Eventplaneinfo *sepds = new Eventplaneinfov1();
      sepds->set_qvector(south_Qvec);
      epmap->insert(sepds, EventplaneinfoMap::sEPDS);

      Eventplaneinfo *sepdn = new Eventplaneinfov1();
      sepdn->set_qvector(north_Qvec);
      epmap->insert(sepdn, EventplaneinfoMap::sEPDN);

      if (Verbosity() > 1) {
        sepds->identify();
        sepdn->identify();
      }
    }

    ResetMe();
  }

  if (_mbdEpReco) {
    ResetMe();

    MbdPmtContainer *mbdpmts =
        findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
    if (!mbdpmts) {
      std::cout << PHWHERE << "::ERROR - cannot find MbdPmtContainer"
                << std::endl;
      exit(-1);
    }

    MbdGeom *mbdgeom = findNode::getClass<MbdGeom>(topNode, "MbdGeom");
    if (!mbdgeom) {
      std::cout << PHWHERE << "::ERROR - cannot find MbdGeom" << std::endl;
      exit(-1);
    }

    if (mbdpmts) {
      if (Verbosity()) {
        std::cout << "EventPlaneReco::process_event -  mbdpmts" << std::endl;
      }

      mbdQ = 0.;

      for (int ipmt = 0; ipmt < mbdpmts->get_npmt(); ipmt++) {
        float mbd_q = mbdpmts->get_pmt(ipmt)->get_q();
        mbdQ += mbd_q;
      }

      for (int ipmt = 0; ipmt < mbdpmts->get_npmt(); ipmt++) {
        float mbd_q = mbdpmts->get_pmt(ipmt)->get_q();
        float phi = mbdgeom->get_phi(ipmt);
        int arm = mbdgeom->get_arm(ipmt);

        if (mbdQ < _mbd_e) {
          continue;
        }

        if (arm == 0) {
          for (unsigned int order = 0; order < m_MaxOrder; order++) {
            double Cosine = cos(phi * (double)(order + 1));
            double Sine = sin(phi * (double)(order + 1));
            south_q[order][0] += mbd_q * Cosine; // south Qn,x
            south_q[order][1] += mbd_q * Sine;   // south Qn,y
          }
        } else if (arm == 1) {
          for (unsigned int order = 0; order < m_MaxOrder; order++) {
            double Cosine = cos(phi * (double)(order + 1));
            double Sine = sin(phi * (double)(order + 1));
            north_q[order][0] += mbd_q * Cosine; // north Qn,x
            north_q[order][1] += mbd_q * Sine;   // north Qn,y
          }
        }
      }
    }

    for (unsigned int order = 0; order < m_MaxOrder; order++) {
      south_Qvec.emplace_back(south_q[order][0], south_q[order][1]);
      north_Qvec.emplace_back(north_q[order][0], north_q[order][1]);
    }

    if (mbdpmts) {
      Eventplaneinfo *mbds = new Eventplaneinfov1();
      mbds->set_qvector(south_Qvec);
      epmap->insert(mbds, EventplaneinfoMap::MBDS);

      Eventplaneinfo *mbdn = new Eventplaneinfov1();
      mbdn->set_qvector(north_Qvec);
      epmap->insert(mbdn, EventplaneinfoMap::MBDN);

      if (Verbosity() > 1) {
        mbds->identify();
        mbdn->identify();
      }
    }

    ResetMe();
  }

  if (Verbosity()) {
    epmap->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventPlaneReco::CreateNodes(PHCompositeNode *topNode) {
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode =
      dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(
      iter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode) {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  EventplaneinfoMap *eps =
      findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!eps) {
    eps = new EventplaneinfoMapv1();
    PHIODataNode<PHObject> *EpMapNode =
        new PHIODataNode<PHObject>(eps, "EventplaneinfoMap", "PHObject");
    globalNode->addNode(EpMapNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void EventPlaneReco::ResetMe() {
  for (auto &vec : south_q) {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  for (auto &vec : north_q) {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  south_Qvec.clear();
  north_Qvec.clear();
}

int EventPlaneReco::End(PHCompositeNode * /*topNode*/) {
  std::cout << " EventPlaneReco::End() " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
