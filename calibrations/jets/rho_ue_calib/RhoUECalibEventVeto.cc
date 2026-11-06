#include "RhoUECalibEventVeto.h"

#include <calotrigger/MinimumBiasInfo.h>

#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdDefs.h>
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

RhoUECalibEventVeto::RhoUECalibEventVeto(const std::string &name)
  : SubsysReco(name)
{
}

int RhoUECalibEventVeto::Init(PHCompositeNode * /*topNode*/)
{
  if (m_main_bunch && !m_t0_table.empty())
  {
    // knee table of RhoUECalibPayload::FitMbdKnees: a header "# run knee ... t0 ..." names the columns
    std::ifstream in(m_t0_table);
    if (!in)
    {
      std::cout << PHWHERE << " " << Name() << ": cannot read t0 table " << m_t0_table << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    int t0_column = -1;
    std::string line;
    while (std::getline(in, line))
    {
      if (line.rfind("# run", 0) == 0)
      {
        std::istringstream ss(line.substr(1));
        std::string col;
        for (int i = 0; ss >> col; ++i)
        {
          if (col == "t0")
          {
            t0_column = i;
          }
        }
        continue;
      }
      if (line.empty() || line[0] == '#' || t0_column < 0)
      {
        continue;
      }
      std::istringstream ss(line);
      std::vector<std::string> cols;
      std::string col;
      while (ss >> col)
      {
        cols.push_back(col);
      }
      if (static_cast<int>(cols.size()) > t0_column)
      {
        const double t0 = std::stod(cols.at(t0_column));
        if (std::isfinite(t0))
        {
          m_t0_offsets[std::stoi(cols.at(0))] = t0;
        }
      }
    }
    if (t0_column < 0)
    {
      std::cout << PHWHERE << " " << Name() << ": no t0 column in " << m_t0_table << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RhoUECalibEventVeto::InitRun(PHCompositeNode * /*topNode*/)
{
  const int run = recoConsts::instance()->get_IntFlag("RUNNUMBER");
  m_t0_offset = m_t0_offset_user;
  m_t0_offset_from_table = false;
  if (m_main_bunch && !m_t0_table.empty())
  {
    const auto it = m_t0_offsets.find(run);
    if (it == m_t0_offsets.end())
    {
      std::cout << PHWHERE << " " << Name() << ": no main-bunch t0 for run " << run << " in " << m_t0_table << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    m_t0_offset = it->second;
    m_t0_offset_from_table = true;
  }
  if (m_main_bunch)
  {
    std::cout << Name() << "::InitRun - run " << run << ": main bunch |t0 - " << m_t0_offset << " ns| < " << m_t0_window
              << " ns" << (m_t0_offset_from_table ? " (table)" : "") << ", or the classifier's bunch class" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RhoUECalibEventVeto::process_event(PHCompositeNode *topNode)
{
  ++m_n_seen;

  auto *mbinfo = findNode::getClass<MinimumBiasInfo>(topNode, m_minbias_node);
  if (m_require_minbias)
  {
    if (!mbinfo)
    {
      std::cout << PHWHERE << " " << Name() << ": MinimumBiasInfo node " << m_minbias_node << " not found" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    if (!mbinfo->isAuAuMinimumBias())
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  ++m_n_minbias;

  // vertex
  if (m_zvtx_max >= 0)
  {
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
    if (!vtx || !std::isfinite(vtx->get_z()) || std::fabs(vtx->get_z()) >= m_zvtx_max)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  ++m_n_vertex;

  auto *mbdout = findNode::getClass<MbdOut>(topNode, m_mbd_node);
  if ((m_main_bunch || m_offband_cut > 0) && !mbdout)
  {
    std::cout << PHWHERE << " " << Name() << ": MBD node " << m_mbd_node << " not found" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // bunch crossing: the classifier's class, or the t0 window around the run's main-bunch centre
  const bool classifier_classes = (mbinfo && std::isfinite(mbinfo->getT0Corrected()));
  if (m_tag_bunch_class && mbinfo && !classifier_classes && mbdout)
  {
    const double t0 = mbdout->get_t0();
    if (std::isfinite(t0))
    {
      const double t0c = t0 - m_t0_offset;
      int bunch = MinimumBiasInfo::BUNCH::NONE;
      if (std::fabs(t0c) < m_main_window)
      {
        bunch = MinimumBiasInfo::BUNCH::MAIN;
      }
      else if (std::fabs(t0c) < m_satellite_window)
      {
        bunch = (t0c < 0) ? MinimumBiasInfo::BUNCH::M5 : MinimumBiasInfo::BUNCH::P5;
      }
      mbinfo->setT0Corrected(static_cast<float>(t0c));
      mbinfo->setBunchClass(bunch);
    }
  }
  if (m_main_bunch)
  {
    bool accept = false;
    if (classifier_classes)
    {
      const int bunch = mbinfo->getBunchClass();
      accept = (bunch == MinimumBiasInfo::BUNCH::MAIN) ||
               (m_accept_satellites && (bunch == MinimumBiasInfo::BUNCH::M5 || bunch == MinimumBiasInfo::BUNCH::P5));
      ++m_n_bunch_classifier;
    }
    else
    {
      const double t0 = mbdout->get_t0();
      accept = std::isfinite(t0) && std::fabs(t0 - m_t0_offset) < m_t0_window;
    }
    if (!accept)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  ++m_n_bunch;

  // off-band: charge lost in the production reconstruction
  if (m_offband_cut > 0)
  {
    auto *pmts = findNode::getClass<MbdPmtContainer>(topNode, m_mbd_pmt_node);
    if (!pmts)
    {
      std::cout << PHWHERE << " " << Name() << ": MBD PMT node " << m_mbd_pmt_node << " not found" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    double qfull = 0;
    for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ++ipmt)
    {
      const float q = pmts->get_pmt(ipmt)->get_q();
      if (q > 0)
      {
        qfull += q;
      }
    }
    MbdOut *prodout = mbdout;
    if (!m_prod_topnode.empty())
    {
      prodout = findNode::getClass<MbdOut>(Fun4AllServer::instance()->topNode(m_prod_topnode), m_mbd_node);
      if (!prodout)
      {
        std::cout << PHWHERE << " " << Name() << ": production MBD node " << m_mbd_node << " not found under " << m_prod_topnode << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }
    }
    const double qprod = prodout->get_q(0) + prodout->get_q(1);
    if (!(qprod > 0))
    {
      ++m_n_no_prod_charge;
    }
    if (qprod > 0 && qfull / qprod > m_offband_cut)
    {
      ++m_n_offband_veto;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  ++m_n_selected;
  return Fun4AllReturnCodes::EVENT_OK;
}

int RhoUECalibEventVeto::End(PHCompositeNode * /*topNode*/)
{
  std::string offband = " (cut off)";
  if (m_offband_cut > 0)
  {
    offband = m_prod_topnode.empty() ? " (vs the refit MbdOut)" : " (vs the production MbdOut, " + std::to_string(m_n_no_prod_charge) + " without production charge)";
  }
  std::cout << Name() << "::End - " << m_n_seen << " events, " << m_n_minbias << " minimum bias" << (m_require_minbias ? "" : " (not required)") << ", " << m_n_vertex << " with a vertex" << (m_zvtx_max >= 0 ? " |z| < " + std::to_string(m_zvtx_max) : "")
            << ", " << m_n_bunch << " in the bunch window" << (m_main_bunch ? " (|t0'| < " + std::to_string(m_t0_window) + " ns, " + std::to_string(m_n_bunch_classifier) + " by the classifier's bunch class)" : " (cut off)")
            << ", " << m_n_offband_veto << " off-band vetoed" << offband << ", " << m_n_selected << " kept" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
