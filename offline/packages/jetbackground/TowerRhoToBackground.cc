#include "TowerRhoToBackground.h"

#include "TowerBackground.h"
#include "TowerBackgroundv1.h"
#include "TowerRho.h"

#include <cdbobjects/CDBTTree.h>

#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdOut.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

TowerRhoToBackground::~TowerRhoToBackground()
{
  delete m_eta_calib;
}

int TowerRhoToBackground::InitRun(PHCompositeNode *topNode)
{
  if (LoadEtaCalib() != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator iter(topNode);
  auto *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << " DST node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  auto *bkgNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "JETBACKGROUND"));
  if (!bkgNode)
  {
    bkgNode = new PHCompositeNode("JETBACKGROUND");
    dstNode->addNode(bkgNode);
  }
  if (!findNode::getClass<TowerBackground>(topNode, m_background_node))
  {
    bkgNode->addNode(new PHIODataNode<PHObject>(new TowerBackgroundv1(), m_background_node, "PHObject"));
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TowerRhoToBackground::process_event(PHCompositeNode *topNode)
{
  auto *background = findNode::getClass<TowerBackground>(topNode, m_background_node);
  if (!background)
  {
    std::cout << PHWHERE << " TowerBackground node " << m_background_node << " not found, exiting" << std::endl;
    exit(1);
  }

  // event vertex z, 0 if there is none
  float vtxz = 0;
  auto *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (vertexmap && !vertexmap->empty())
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if (m_use_vertex_type)
    {
      const auto vertices = vertexmap->get_gvtxs_with_type({m_vertex_type});
      vtx = vertices.empty() ? nullptr : vertices.at(0);
    }
    if (vtx && std::isfinite(vtx->get_z()))
    {
      vtxz = vtx->get_z();
    }
  }

  ++m_n_events;

  // calibration bin of this event, -1 (w = 1) if no calibration or out of range
  int izbin = -1;
  int imbd = -1;
  if (m_eta_calib)
  {
    auto *mbdout = findNode::getClass<MbdOut>(topNode, m_mbd_node);
    if (!mbdout)
    {
      std::cout << PHWHERE << " MBD node " << m_mbd_node << " needed for the eta calibration not found, exiting" << std::endl;
      exit(1);
    }
    const float mbdq = mbdout->get_q(0) + mbdout->get_q(1);
    izbin = find_bin(vtxz, m_calib_zvtx_edges);
    imbd = find_bin(mbdq, m_calib_mbdq_edges);
    if (izbin < 0 || imbd < 0)
    {
      ++m_n_out_of_range;
      m_n_zvtx_out += (izbin < 0) ? 1 : 0;
      m_n_mbdq_out += (imbd < 0) ? 1 : 0;
      if (Verbosity() > 1)
      {
        std::cout << "TowerRhoToBackground::process_event - fallback to w = 1: vertex z = " << vtxz
                  << (izbin < 0 ? " (out of range)" : "") << ", MBD charge = " << mbdq
                  << (imbd < 0 ? " (out of range)" : "") << std::endl;
      }
    }
  }
  else
  {
    ++m_n_no_calib;
  }
  bool event_has_invalid_w = false;
  const int n_mbd_bins = static_cast<int>(m_calib_mbdq_edges.size()) - 1;

  // layer 0 is the retowered EMCal: IHCal strips at the EMCal radius
  const RawTowerDefs::CalorimeterId caloid[3] = {RawTowerDefs::CEMC, RawTowerDefs::HCALIN, RawTowerDefs::HCALOUT};
  const int strip_layer[3] = {1, 1, 2};
  const std::string calib_field[3] = {"w_cemc", "w_hcalin", "w_hcalout"};

  for (int layer = 0; layer < 3; layer++)
  {
    auto *towerrho = findNode::getClass<TowerRho>(topNode, m_rho_nodes[layer]);
    auto *radius_geom = findNode::getClass<RawTowerGeomContainer>(topNode, m_geom_nodes[layer]);
    auto *strip_geom = findNode::getClass<RawTowerGeomContainer>(topNode, m_geom_nodes[strip_layer[layer]]);
    if (!towerrho || !radius_geom || !strip_geom)
    {
      std::cout << PHWHERE << " missing " << m_rho_nodes[layer] << " or tower geometry, exiting" << std::endl;
      exit(1);
    }

    const double radius = radius_geom->get_tower_geometry(RawTowerDefs::encode_towerid(caloid[layer], 0, 0))->get_center_radius();
    auto corrected_eta = [radius, vtxz](const double eta)
    { return std::asinh(((std::sinh(eta) * radius) - vtxz) / radius); };

    const bool is_area = (towerrho->get_method() == TowerRho::Method::AREA);
    const auto phibounds = strip_geom->get_phibounds(0);
    const double dphi = phibounds.second - phibounds.first;

    if (m_eta_calib && m_calib_neta != strip_geom->get_etabins())
    {
      std::cout << PHWHERE << " eta calibration has " << m_calib_neta << " eta bins, "
                << m_geom_nodes[strip_layer[layer]] << " has " << strip_geom->get_etabins() << ", exiting" << std::endl;
      exit(1);
    }

    std::vector<float> ue(strip_geom->get_etabins(), 0);
    for (int ieta = 0; ieta < strip_geom->get_etabins(); ieta++)
    {
      double ue_tower = towerrho->get_rho() * std::cosh(corrected_eta(strip_geom->get_etacenter(ieta)));
      if (is_area)
      {
        const auto etabounds = strip_geom->get_etabounds(ieta);
        ue_tower *= (corrected_eta(etabounds.second) - corrected_eta(etabounds.first)) * dphi;
      }
      if (izbin >= 0 && imbd >= 0)
      {
        const float w = m_eta_calib->GetFloatValue((izbin * n_mbd_bins * m_calib_neta) + (imbd * m_calib_neta) + ieta, calib_field[layer]);
        if (w > 0)
        {
          ue_tower *= w;
        }
        else
        {
          ++m_n_invalid_w;
          event_has_invalid_w = true;
          if (Verbosity() > 1)
          {
            std::cout << "TowerRhoToBackground::process_event - fallback to w = 1: invalid weight " << w
                      << " for " << calib_field[layer] << ", eta strip " << ieta << " (vertex z bin " << izbin
                      << ", MBD charge bin " << imbd << ")" << std::endl;
          }
        }
      }
      ue[ieta] = static_cast<float>(ue_tower);
    }
    background->set_UE(layer, ue);

    if (Verbosity() > 1)
    {
      std::cout << "TowerRhoToBackground::process_event - " << m_rho_nodes[layer] << " rho = " << towerrho->get_rho()
                << (is_area ? " (AREA)" : " (MULT)") << ", vertex z = " << vtxz << std::endl;
    }
  }
  m_n_events_invalid_w += event_has_invalid_w ? 1 : 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

int TowerRhoToBackground::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "TowerRhoToBackground::End - " << m_n_events << " events";
  if (m_n_no_calib > 0)
  {
    std::cout << "; no eta calibration configured, w = 1 in all " << m_n_no_calib << " of them";
  }
  std::cout << std::endl;
  if (m_n_events > m_n_no_calib)
  {
    std::cout << "TowerRhoToBackground::End - fallbacks to w = 1: " << m_n_out_of_range
              << " events out of the calibrated range (vertex z: " << m_n_zvtx_out
              << ", MBD charge: " << m_n_mbdq_out << "), " << m_n_invalid_w << " invalid weights in "
              << m_n_events_invalid_w << " events" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TowerRhoToBackground::LoadEtaCalib()
{
  delete m_eta_calib;
  m_eta_calib = nullptr;
  if (m_eta_calib_path.empty() && m_eta_calib_tag.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;  // no calibration, w = 1
  }

  const std::string url = m_eta_calib_path.empty() ? CDBInterface::instance()->getUrl(m_eta_calib_tag) : m_eta_calib_path;
  if (url.empty())
  {
    std::cout << PHWHERE << " no eta calibration found for CDB tag " << m_eta_calib_tag << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_eta_calib = new CDBTTree(url);
  m_eta_calib->LoadCalibrations();
  m_calib_neta = m_eta_calib->GetSingleIntValue("n_eta");
  const int n_zvtx = m_eta_calib->GetSingleIntValue("n_zvtx_bins");
  const int n_mbdq = m_eta_calib->GetSingleIntValue("n_mbdQ_bins");
  if (m_calib_neta <= 0 || n_zvtx <= 0 || n_mbdq <= 0)
  {
    std::cout << PHWHERE << " invalid eta calibration " << url << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_calib_zvtx_edges.assign(n_zvtx + 1, 0);
  for (int i = 0; i <= n_zvtx; i++)
  {
    m_calib_zvtx_edges[i] = m_eta_calib->GetSingleFloatValue("zvtx_edge_" + std::to_string(i));
  }
  m_calib_mbdq_edges.assign(n_mbdq + 1, 0);
  for (int i = 0; i <= n_mbdq; i++)
  {
    m_calib_mbdq_edges[i] = m_eta_calib->GetSingleFloatValue("mbdQ_edge_" + std::to_string(i));
  }

  // bin edges must be finite and strictly increasing
  for (const auto *edges : {&m_calib_zvtx_edges, &m_calib_mbdq_edges})
  {
    for (size_t i = 0; i < edges->size(); i++)
    {
      if (!std::isfinite(edges->at(i)) || (i > 0 && edges->at(i) <= edges->at(i - 1)))
      {
        std::cout << PHWHERE << " invalid bin edges in eta calibration " << url << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << "TowerRhoToBackground::LoadEtaCalib - " << url << " (" << m_calib_neta << " eta x "
              << n_zvtx << " vertex z x " << n_mbdq << " MBD charge bins)" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TowerRhoToBackground::find_bin(const float val, const std::vector<float> &edges)
{
  for (size_t i = 0; i + 1 < edges.size(); i++)
  {
    if (val >= edges[i] && val < edges[i + 1])
    {
      return static_cast<int>(i);
    }
  }
  return -1;
}
