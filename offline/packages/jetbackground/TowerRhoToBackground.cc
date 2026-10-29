#include "TowerRhoToBackground.h"

#include "MbdChargeKnee.h"
#include "TowerBackground.h"
#include "TowerBackgroundv1.h"
#include "TowerRho.h"

#include <cdbobjects/CDBHistos.h>

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
#include <phool/recoConsts.h>

#include <TH1.h>
#include <TH3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

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
    std::cout << PHWHERE << " TowerBackground node " << m_background_node << " not found" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // event vertex z, 0 if there is none
  float vtxz = 0;
  bool have_vertex = false;
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
      have_vertex = true;
    }
  }

  ++m_n_events;

  // calibration nodes of this event in (z, Q'); no calibration (UE = x) if none is
  // configured, without a vertex, or outside the calibrated vertex range
  bool calibrate = m_calib_on;
  Nodes znodes;
  Nodes qnodes;
  if (!m_calib_on)
  {
    ++m_n_no_calib;
  }
  else if (!have_vertex || vtxz < m_zvtx_edges.front() || vtxz >= m_zvtx_edges.back())
  {
    calibrate = false;
    ++m_n_out_of_range;
    if (Verbosity() > 1)
    {
      std::cout << "TowerRhoToBackground::process_event - no calibration: "
                << (have_vertex ? "vertex z = " + std::to_string(vtxz) + " out of range" : std::string("no vertex")) << std::endl;
    }
  }
  else
  {
    auto *mbdout = findNode::getClass<MbdOut>(topNode, m_mbd_node);
    if (!mbdout)
    {
      std::cout << PHWHERE << " MBD node " << m_mbd_node << " needed for the calibration not found" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    const float qprime = (mbdout->get_q(0) + mbdout->get_q(1)) / m_knee;
    m_n_zvtx_clamped += find_nodes(vtxz, m_zvtx_centres, m_zvtx_edges, m_interpolate, znodes) ? 1 : 0;
    m_n_mbdq_clamped += find_nodes(qprime, m_mbdq_centres, m_mbdq_edges, m_interpolate, qnodes) ? 1 : 0;
  }
  bool event_has_invalid_w = false;
  const int n_mbdq = static_cast<int>(m_mbdq_centres.size());

  // bilinear value of a calibration table at this event's nodes
  auto value = [&](const std::vector<float> &table, const int ieta)
  {
    auto at = [&](const int iz, const int iq)
    { return table[(((iz * n_mbdq) + iq) * m_calib_neta) + ieta]; };
    return ((1 - znodes.f) * (1 - qnodes.f) * at(znodes.i0, qnodes.i0)) + ((1 - znodes.f) * qnodes.f * at(znodes.i0, qnodes.i1)) +
           (znodes.f * (1 - qnodes.f) * at(znodes.i1, qnodes.i0)) + (znodes.f * qnodes.f * at(znodes.i1, qnodes.i1));
  };

  // layer 0 is the retowered EMCal: IHCal strips at the EMCal radius
  const RawTowerDefs::CalorimeterId caloid[3] = {RawTowerDefs::CEMC, RawTowerDefs::HCALIN, RawTowerDefs::HCALOUT};
  const int strip_layer[3] = {1, 1, 2};

  for (int layer = 0; layer < 3; layer++)
  {
    auto *towerrho = findNode::getClass<TowerRho>(topNode, m_rho_nodes[layer]);
    auto *radius_geom = findNode::getClass<RawTowerGeomContainer>(topNode, m_geom_nodes[layer]);
    auto *strip_geom = findNode::getClass<RawTowerGeomContainer>(topNode, m_geom_nodes[strip_layer[layer]]);
    if (!towerrho || !radius_geom || !strip_geom)
    {
      std::cout << PHWHERE << " missing " << m_rho_nodes[layer] << " or tower geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    const double radius = radius_geom->get_tower_geometry(RawTowerDefs::encode_towerid(caloid[layer], 0, 0))->get_center_radius();
    auto corrected_eta = [radius, vtxz](const double eta)
    { return std::asinh(((std::sinh(eta) * radius) - vtxz) / radius); };

    const bool is_area = (towerrho->get_method() == TowerRho::Method::AREA);
    const auto phibounds = strip_geom->get_phibounds(0);
    const double dphi = phibounds.second - phibounds.first;

    if (m_calib_on && m_calib_neta != strip_geom->get_etabins())
    {
      std::cout << PHWHERE << " calibration has " << m_calib_neta << " eta bins, "
                << m_geom_nodes[strip_layer[layer]] << " has " << strip_geom->get_etabins() << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
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
      if (calibrate)
      {
        bool valid = false;
        if (m_model == CalibModel::AFFINE)
        {
          const float alpha = value(m_alpha[layer], ieta);
          const float beta = value(m_beta[layer], ieta);
          valid = std::isfinite(alpha) && std::isfinite(beta);
          if (valid)
          {
            ue_tower = (alpha * ue_tower) + beta;
          }
        }
        else
        {
          const float w = value(m_w[layer], ieta);
          valid = std::isfinite(w) && w > 0;
          if (valid)
          {
            ue_tower *= w;
          }
        }
        if (!valid)
        {
          ++m_n_invalid_w;
          event_has_invalid_w = true;
          if (Verbosity() > 1)
          {
            std::cout << "TowerRhoToBackground::process_event - invalid calibration for layer " << layer
                      << ", eta strip " << ieta << ": UE = rho x cosh(eta)" << std::endl;
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
    std::cout << "; no calibration configured, UE = rho x cosh(eta) in all " << m_n_no_calib << " of them";
  }
  std::cout << std::endl;
  if (m_n_events > m_n_no_calib)
  {
    std::cout << "TowerRhoToBackground::End - " << (m_model == CalibModel::AFFINE ? "affine" : "multiplicative")
              << " calibration; no vertex or vertex out of range (uncalibrated): " << m_n_out_of_range
              << "; clamped beyond the outermost bin centres: vertex z " << m_n_zvtx_clamped << ", Q' " << m_n_mbdq_clamped
              << "; " << m_n_invalid_w << " invalid calibrations in " << m_n_events_invalid_w << " events" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TowerRhoToBackground::LoadEtaCalib()
{
  m_calib_on = false;
  for (auto *table : {&m_w, &m_alpha, &m_beta})
  {
    for (auto &layer : *table)
    {
      layer.clear();
    }
  }
  if (m_eta_calib_path.empty() && m_eta_calib_tag.empty())
  {
    return Fun4AllReturnCodes::EVENT_OK;  // no calibration, UE = x
  }

  const std::string url = m_eta_calib_path.empty() ? CDBInterface::instance()->getUrl(m_eta_calib_tag) : m_eta_calib_path;
  if (url.empty())
  {
    std::cout << PHWHERE << " no calibration found for CDB tag " << m_eta_calib_tag << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // the payload is only needed here: copy it into flat vectors and release it
  CDBHistos cdb(url);
  cdb.LoadCalibrations();
  TH1 *hzvtx = cdb.getHisto("zvtx_bins");
  TH1 *hmbdq = cdb.getHisto("mbdq_bins");
  TH1 *hknee = cdb.getHisto("mbdq_knee", false);
  if (!hzvtx || !hmbdq)
  {
    std::cout << PHWHERE << " calibration " << url << ": zvtx_bins or mbdq_bins missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // bin edges (finite, strictly increasing) and centres (inside their bins)
  const std::array<std::pair<TH1 *, std::pair<std::vector<float> *, std::vector<float> *>>, 2> axes = {
    std::make_pair(hzvtx, std::make_pair(&m_zvtx_edges, &m_zvtx_centres)),
    std::make_pair(hmbdq, std::make_pair(&m_mbdq_edges, &m_mbdq_centres))};
  for (const auto &[h, vecs] : axes)
  {
    const int nbins = h->GetNbinsX();
    auto &edges = *vecs.first;
    auto &centres = *vecs.second;
    edges.assign(nbins + 1, 0);
    centres.assign(nbins, 0);
    for (int i = 0; i <= nbins; i++)
    {
      edges[i] = static_cast<float>(h->GetXaxis()->GetBinLowEdge(i + 1));
    }
    for (int i = 0; i < nbins; i++)
    {
      centres[i] = static_cast<float>(h->GetBinContent(i + 1));
    }
    for (int i = 0; i <= nbins; i++)
    {
      const bool bad_edge = !std::isfinite(edges[i]) || (i > 0 && edges[i] <= edges[i - 1]);
      const bool bad_centre = (i < nbins) && !(centres[i] >= edges[i] && centres[i] <= h->GetXaxis()->GetBinUpEdge(i + 1));
      if (bad_edge || bad_centre)
      {
        std::cout << PHWHERE << " invalid bin edges or centres in calibration " << url << " (" << h->GetName() << ")" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
  }

  // Q_knee of this run
  const int run = recoConsts::instance()->get_IntFlag("RUNNUMBER");
  std::string knee_source;
  bool knee_default = false;
  m_knee = MbdChargeKnee::find(m_knee_user, m_knee_tag, m_knee_path, hknee, run, knee_source, knee_default);
  if (!(m_knee > 0))
  {
    std::cout << PHWHERE << " UE calibration: no MBD charge knee for run " << run << ": " << knee_source << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (knee_default)
  {
    std::cout << "TowerRhoToBackground::LoadEtaCalib - WARNING run " << run << " Q_knee = " << m_knee << " (" << knee_source << ")" << std::endl;
  }

  // tables per layer: TH3 (ieta, z, Q')
  const int n_zvtx = static_cast<int>(m_zvtx_centres.size());
  const int n_mbdq = static_cast<int>(m_mbdq_centres.size());
  const std::array<std::string, 3> layer_name = {"CEMC", "HCALIN", "HCALOUT"};
  m_calib_neta = 0;
  for (int layer = 0; layer < 3; layer++)
  {
    const std::array<std::pair<std::string, std::vector<float> *>, 3> tables = {
      std::make_pair("w_" + layer_name[layer], &m_w[layer]),
      std::make_pair("alpha_" + layer_name[layer], &m_alpha[layer]),
      std::make_pair("beta_" + layer_name[layer], &m_beta[layer])};
    for (const auto &[name, table] : tables)
    {
      auto *h = dynamic_cast<TH3 *>(cdb.getHisto(name));
      if (!h || h->GetNbinsY() != n_zvtx || h->GetNbinsZ() != n_mbdq || (m_calib_neta > 0 && h->GetNbinsX() != m_calib_neta))
      {
        std::cout << PHWHERE << " calibration " << url << ": " << name << " missing or not n_eta x "
                  << n_zvtx << " x " << n_mbdq << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
      m_calib_neta = h->GetNbinsX();
      table->assign(static_cast<size_t>(n_zvtx) * n_mbdq * m_calib_neta, 0);
      for (int iz = 0; iz < n_zvtx; iz++)
      {
        for (int iq = 0; iq < n_mbdq; iq++)
        {
          for (int ieta = 0; ieta < m_calib_neta; ieta++)
          {
            (*table)[(((iz * n_mbdq) + iq) * m_calib_neta) + ieta] = static_cast<float>(h->GetBinContent(ieta + 1, iz + 1, iq + 1));
          }
        }
      }
    }
  }
  m_calib_on = true;

  if (Verbosity() > 0)
  {
    std::cout << "TowerRhoToBackground::LoadEtaCalib - " << url << " (" << m_calib_neta << " eta x " << n_zvtx
              << " vertex z x " << n_mbdq << " Q' bins), run " << run << " Q_knee = " << m_knee << " (" << knee_source << "), "
              << (m_model == CalibModel::AFFINE ? "affine" : "multiplicative")
              << (m_interpolate ? ", interpolated" : ", bin lookup") << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

bool TowerRhoToBackground::find_nodes(const float val, const std::vector<float> &centres, const std::vector<float> &edges,
                                      const bool interpolate, Nodes &nodes)
{
  const int nbins = static_cast<int>(centres.size());
  if (!interpolate)
  {
    // the bin of val, the first / last bin beyond the edges
    const auto it = std::upper_bound(edges.begin(), edges.end(), val);
    const int ibin = std::clamp(static_cast<int>(it - edges.begin()) - 1, 0, nbins - 1);
    nodes = {ibin, ibin, 0};
    return val < edges.front() || val >= edges.back();
  }
  if (val <= centres.front())
  {
    nodes = {0, 0, 0};
    return val < centres.front();
  }
  if (val >= centres.back())
  {
    nodes = {nbins - 1, nbins - 1, 0};
    return val > centres.back();
  }
  const auto it = std::upper_bound(centres.begin(), centres.end(), val);
  const int i1 = static_cast<int>(it - centres.begin());
  const int i0 = i1 - 1;
  nodes = {i0, i1, (val - centres[i0]) / (centres[i1] - centres[i0])};
  return false;
}
