#include "RhoUECalibBase.h"

#include <calotrigger/MinimumBiasInfo.h>

#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdOut.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

RhoUECalibBase::RhoUECalibBase(const std::string &name, const std::string &outfile)
  : SubsysReco(name)
  , m_outfile(outfile)
{
}

RhoUECalibBase::~RhoUECalibBase()
{
  delete m_file;
}

std::vector<double> RhoUECalibBase::default_zvtx_edges()
{
  // 2.5 cm to |z| = 45 cm, then 45-52.5 and 52.5-60 (cross-validated, run 2024 Au+Au)
  std::vector<double> edges = {-60, -52.5};
  for (int i = 0; i <= 36; ++i)
  {
    edges.push_back(-45 + (2.5 * i));
  }
  edges.push_back(52.5);
  edges.push_back(60);
  return edges;
}

std::vector<double> RhoUECalibBase::default_mbdq_edges()
{
  // fine at low Q' (cross-validated, run 2024 Au+Au)
  return {0, 0.01, 0.03, 0.055, 0.085, 0.115, 0.15, 0.19, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475,
          0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925, 0.975, 1.025, 1.1, 1.3};
}

const std::array<std::string, RhoUECalibBase::NLAYER> &RhoUECalibBase::layer_names()
{
  static const std::array<std::string, NLAYER> names = {"CEMC", "HCALIN", "HCALOUT"};
  return names;
}

std::string RhoUECalibBase::tower_node(const int layer) const
{
  static const std::array<std::string, NLAYER> suffix = {"_CEMC_RETOWER", "_HCALIN", "_HCALOUT"};
  return m_tower_prefix + suffix.at(layer);
}

int RhoUECalibBase::Init(PHCompositeNode * /*topNode*/)
{
  for (const auto *edges : {&m_zvtx_edges, &m_mbdq_edges})
  {
    if (edges->size() < 2 || !std::is_sorted(edges->begin(), edges->end()) ||
        std::adjacent_find(edges->begin(), edges->end()) != edges->end())
    {
      std::cout << PHWHERE << " " << Name() << ": class bin edges must be strictly increasing" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  if (!m_knee_table.empty())
  {
    std::ifstream in(m_knee_table);
    if (!in)
    {
      std::cout << PHWHERE << " " << Name() << ": cannot read knee table " << m_knee_table << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    std::string line;
    while (std::getline(in, line))
    {
      if (line.empty() || line[0] == '#')
      {
        continue;
      }
      std::istringstream ss(line);
      int run = 0;
      double knee = 0;
      if (ss >> run >> knee && knee > 0)
      {
        m_knees[run] = knee;
      }
    }
  }

  TDirectory *olddir = gDirectory;
  m_file = new TFile(m_outfile.c_str(), "RECREATE");
  auto *hsel = new TH1D("h_sel", "event selection;;events", 4, 0, 4);
  hsel->GetXaxis()->SetBinLabel(1, "seen");
  hsel->GetXaxis()->SetBinLabel(2, "minimum bias");
  hsel->GetXaxis()->SetBinLabel(3, "vertex in range");
  hsel->GetXaxis()->SetBinLabel(4, "selected");
  // class binning (the axes carry the edges; contents unused)
  new TH1D("zvtx_bins", "z_{vtx} class bins;z_{vtx} [cm]", n_zvtx(), m_zvtx_edges.data());
  new TH1D("mbdq_bins", "Q' class bins;Q' = Q/Q_{knee}", n_mbdq(), m_mbdq_edges.data());
  book();
  olddir->cd();
  return Fun4AllReturnCodes::EVENT_OK;
}

int RhoUECalibBase::InitRun(PHCompositeNode * /*topNode*/)
{
  m_run = recoConsts::instance()->get_IntFlag("RUNNUMBER");
  m_knee = m_knee_user;
  if (!(m_knee > 0))
  {
    const auto it = m_knees.find(m_run);
    m_knee = (it != m_knees.end()) ? it->second : -1;
  }
  if (needs_class() && !(m_knee > 0))
  {
    std::cout << PHWHERE << " " << Name() << ": no MBD charge knee for run " << m_run
              << " (set_mbd_knee or set_knee_table)" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // knee used for this run, hadd-safe: bin 1 = sum of the knee over selected events, bin 2 = events
  if (needs_class())
  {
    const std::string name = "knee_" + std::to_string(m_run);
    if (!m_file->Get(name.c_str()))
    {
      TDirectory *olddir = gDirectory;
      m_file->cd();
      new TH1D(name.c_str(), ("MBD charge knee of run " + std::to_string(m_run) + ": knee x events, events").c_str(), 2, 0, 2);
      olddir->cd();
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RhoUECalibBase::process_event(PHCompositeNode *topNode)
{
  ++m_n_seen;
  if (m_require_minbias)
  {
    auto *mbinfo = findNode::getClass<MinimumBiasInfo>(topNode, m_minbias_node);
    if (!mbinfo)
    {
      std::cout << PHWHERE << " " << Name() << ": MinimumBiasInfo node " << m_minbias_node << " not found" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    if (!mbinfo->isAuAuMinimumBias())
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }
  ++m_n_minbias;

  // vertex inside the class range in z
  auto *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  GlobalVertex *vtx = nullptr;
  if (vertexmap && !vertexmap->empty())
  {
    vtx = vertexmap->begin()->second;
    if (m_use_vertex_type)
    {
      const auto vertices = vertexmap->get_gvtxs_with_type({m_vertex_type});
      vtx = vertices.empty() ? nullptr : vertices.at(0);
    }
  }
  if (!vtx || !std::isfinite(vtx->get_z()) || vtx->get_z() < m_zvtx_edges.front() || vtx->get_z() >= m_zvtx_edges.back())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  ++m_n_vertex;
  m_zvtx = vtx->get_z();

  auto *mbdout = findNode::getClass<MbdOut>(topNode, m_mbd_node);
  if (!mbdout)
  {
    std::cout << PHWHERE << " " << Name() << ": MBD node " << m_mbd_node << " not found" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_mbdq = mbdout->get_q(0) + mbdout->get_q(1);
  m_mbdt0 = mbdout->get_t0();

  if (needs_class())
  {
    m_qprime = m_mbdq / m_knee;
    m_izvtx = static_cast<int>(std::upper_bound(m_zvtx_edges.begin(), m_zvtx_edges.end(), m_zvtx) - m_zvtx_edges.begin()) - 1;
    m_imbdq = std::clamp(static_cast<int>(std::upper_bound(m_mbdq_edges.begin(), m_mbdq_edges.end(), m_qprime) - m_mbdq_edges.begin()) - 1,
                         0, n_mbdq() - 1);
    m_class = (m_izvtx * n_mbdq()) + m_imbdq;
    auto *hknee = static_cast<TH1D *>(m_file->Get(("knee_" + std::to_string(m_run)).c_str()));
    hknee->AddBinContent(1, m_knee);
    hknee->AddBinContent(2, 1);
  }
  if (needs_geometry())
  {
    if (!m_geom_loaded && load_geometry(topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
    for (int layer = 0; layer < NLAYER; ++layer)
    {
      const double r = m_radius.at(layer);
      for (int ieta = 0; ieta < NETA; ++ieta)
      {
        m_eta_corr.at(layer).at(ieta) = std::asinh(((std::sinh(m_strip_eta.at(layer).at(ieta)) * r) - m_zvtx) / r);
      }
    }
  }
  ++m_n_selected;
  m_evt = static_cast<int>(m_n_selected);

  TDirectory *olddir = gDirectory;
  m_file->cd();
  const int ret = fill(topNode);
  olddir->cd();
  return ret;
}

int RhoUECalibBase::load_geometry(PHCompositeNode *topNode)
{
  // layer 0 is the retowered EMCal: IHCal strips at the EMCal radius
  const std::array<RawTowerDefs::CalorimeterId, NLAYER> caloid = {RawTowerDefs::CEMC, RawTowerDefs::HCALIN, RawTowerDefs::HCALOUT};
  const std::array<int, NLAYER> strip_layer = {1, 1, 2};
  for (int layer = 0; layer < NLAYER; ++layer)
  {
    auto *radius_geom = findNode::getClass<RawTowerGeomContainer>(topNode, m_geom_nodes.at(layer));
    auto *strip_geom = findNode::getClass<RawTowerGeomContainer>(topNode, m_geom_nodes.at(strip_layer.at(layer)));
    if (!radius_geom || !strip_geom || strip_geom->get_etabins() != NETA)
    {
      std::cout << PHWHERE << " " << Name() << ": tower geometry " << m_geom_nodes.at(layer) << " / "
                << m_geom_nodes.at(strip_layer.at(layer)) << " missing or not " << NETA << " eta bins" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    auto *tower0 = radius_geom->get_tower_geometry(RawTowerDefs::encode_towerid(caloid.at(layer), 0, 0));
    if (!tower0)
    {
      std::cout << PHWHERE << " " << Name() << ": no tower (0, 0) in " << m_geom_nodes.at(layer) << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    m_radius.at(layer) = tower0->get_center_radius();
    for (int ieta = 0; ieta < NETA; ++ieta)
    {
      m_strip_eta.at(layer).at(ieta) = strip_geom->get_etacenter(ieta);
    }
  }
  m_geom_loaded = true;
  return Fun4AllReturnCodes::EVENT_OK;
}

int RhoUECalibBase::End(PHCompositeNode * /*topNode*/)
{
  auto *hsel = static_cast<TH1D *>(m_file->Get("h_sel"));
  hsel->SetBinContent(1, static_cast<double>(m_n_seen));
  hsel->SetBinContent(2, static_cast<double>(m_n_minbias));
  hsel->SetBinContent(3, static_cast<double>(m_n_vertex));
  hsel->SetBinContent(4, static_cast<double>(m_n_selected));
  m_file->Write();
  m_file->Close();
  std::cout << Name() << "::End - " << m_n_seen << " events, " << m_n_minbias << " minimum bias, " << m_n_vertex
            << " with a vertex in range, " << m_n_selected << " selected -> " << m_outfile << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
