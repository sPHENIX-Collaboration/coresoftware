#include "CentralityReco.h"

#include "CentralityInfov2.h"

#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <calotrigger/MinimumBiasInfo.h>

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

namespace
{
  // class-suffixed payload columns, indexed by MinimumBiasInfo::BUNCH
  const std::array<std::string, 3> bunch_suffix{"", "_m5", "_p5"};
}  // namespace

CentralityReco::CentralityReco(const std::string &name)
  : SubsysReco(name)
{
}

int CentralityReco::InitRun(PHCompositeNode *topNode)
{
  CDBInterface *m_cdb = CDBInterface::instance();

  // the scale payload decides whether this run uses bunch classes, so it is read first
  std::string centscale_url;
  if (!m_overwrite_url_scale.empty())
  {
    centscale_url = m_overwrite_url_scale;
    std::cout << " Overwriting Scale to " << m_overwrite_url_scale << std::endl;
  }
  else
  {
    centscale_url = m_cdb->getUrl("CentralityScale");
  }

  if (Download_centralityScale(centscale_url))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::string centdiv_url;
  if (!m_overwrite_url_divs.empty())
  {
    centdiv_url = m_overwrite_url_divs;
    std::cout << " Overwriting Divs to " << m_overwrite_url_divs << std::endl;
  }
  else
  {
    centdiv_url = m_cdb->getUrl("Centrality");
  }

  if (Download_centralityDivisions(centdiv_url))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::string vertexscale_url;
  if (!m_overwrite_url_vtx.empty())
  {
    vertexscale_url = m_overwrite_url_vtx;
    std::cout << " Overwriting Vtx to " << m_overwrite_url_vtx << std::endl;
  }
  else
  {
    vertexscale_url = m_cdb->getUrl("CentralityVertexScale");
  }

  if (Download_centralityVertexScales(vertexscale_url))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CreateNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::Download_centralityScale(const std::string &dbfile)
{
  m_centrality_scale.fill(1.00);
  m_bunch_classes = false;
  m_t0_offset = 0.;

  std::filesystem::path dbase_file = dbfile;
  if (dbase_file.extension() == ".root")
  {
    CDBTTree *cdbttree = new CDBTTree(dbase_file);
    cdbttree->LoadCalibrations();
    m_centrality_scale[0] = cdbttree->GetDoubleValue(0, "centralityscale");
    // bunch-class payload: t0 offset, bunch spacings and per-class scales in the same file
    if (m_use_bunch_classes && !m_sum_all_channels)
    {
      const double t0_offset = cdbttree->GetDoubleValue(0, "t0offset");
      if (std::isfinite(t0_offset))
      {
        m_bunch_classes = true;
        m_t0_offset = t0_offset;
        for (int ib = 1; ib < 3; ib++)
        {
          m_centrality_scale[ib] = cdbttree->GetDoubleValue(0, "centralityscale" + bunch_suffix[ib]);
        }
      }
    }
    if (Verbosity())
    {
      std::cout << "centscale = " << m_centrality_scale[0] << " " << m_centrality_scale[1] << " " << m_centrality_scale[2]
                << (m_bunch_classes ? " (bunch classes)" : " (single class)") << " t0offset = " << m_t0_offset << std::endl;
    }
    delete cdbttree;
  }
  else
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbfile << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::Download_centralityDivisions(const std::string &dbfile)
{
  for (auto &map : m_centrality_map)
  {
    map.fill(0);
  }

  std::filesystem::path dbase_file = dbfile;

  if (dbase_file.extension() == ".root")
  {
    CDBTTree *cdbttree = new CDBTTree(dbase_file);
    cdbttree->LoadCalibrations();
    if (Verbosity())
    {
      cdbttree->Print();
    }
    const int nclass = m_bunch_classes ? 3 : 1;
    for (int ib = 0; ib < nclass; ib++)
    {
      for (int idiv = 0; idiv < NDIVS; idiv++)
      {
        m_centrality_map[ib][idiv] = cdbttree->GetFloatValue(idiv, "centralitydiv" + bunch_suffix[ib]);
        if (Verbosity())
        {
          std::cout << "centdiv" << bunch_suffix[ib] << " " << idiv << " : " << m_centrality_map[ib][idiv] << std::endl;
        }
      }
    }
    delete cdbttree;
  }
  else
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbfile << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
int CentralityReco::Download_centralityVertexScales(const std::string &dbfile)
{
  std::filesystem::path dbase_file = dbfile;

  if (dbase_file.extension() == ".root")
  {
    CDBTTree *cdbttree = new CDBTTree(dbase_file);
    cdbttree->LoadCalibrations();
    if (Verbosity())
    {
      cdbttree->Print();
    }

    int nvertexbins = cdbttree->GetIntValue(0, "nvertexbins");

    const int nclass = m_bunch_classes ? 3 : 1;
    for (auto &scales : m_vertex_scales)
    {
      scales.clear();
    }
    for (int ib = 0; ib < nclass; ib++)
    {
      for (int iv = 0; iv < nvertexbins; iv++)
      {
        float scale = cdbttree->GetDoubleValue(iv, "scale" + bunch_suffix[ib]);
        float lowvertex = cdbttree->GetDoubleValue(iv, "low_vertex");
        float highvertex = cdbttree->GetDoubleValue(iv, "high_vertex");
        m_vertex_scales[ib].emplace_back(std::make_pair(lowvertex, highvertex), scale);
      }
    }

    delete cdbttree;
  }
  else
  {
    std::cout << PHWHERE << ", ERROR, unknown file type, " << dbfile << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::ResetEvent(PHCompositeNode * /*unused*/)
{
  m_mbd_total_charge = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::FillVars()
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  const int bunch = m_bunch_classes ? m_mb_info->getBunchClass() : MinimumBiasInfo::BUNCH::MAIN;
  float scale_factor = getVertexScale(bunch);

  if (Verbosity())
  {
    std::cout << scale_factor << "*" << m_centrality_scale[bunch] << std::endl;
  }

  for (int i = 0; i < 128; i++)
  {
    m_mbd_hit = m_mbd_container->get_pmt(i);

    if (m_bunch_classes)
    {
      if (m_mbd_hit->get_q() <= mbd_charge_cut || std::fabs(m_mbd_hit->get_time() - m_t0_offset) >= m_hit_time_window)
      {
        continue;
      }
    }
    else if (!m_sum_all_channels)
    {
      if ((m_mbd_hit->get_q()) < mbd_charge_cut)
      {
        continue;
      }
      if (fabs(m_mbd_hit->get_time()) > mbd_time_cut)
      {
        continue;
      }
    }
    m_mbd_total_charge += m_mbd_hit->get_q() * scale_factor * m_centrality_scale[bunch];
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::FillCentralityInfo()
{
  // Fill is minbias

  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  const int bunch = m_bunch_classes ? m_mb_info->getBunchClass() : MinimumBiasInfo::BUNCH::MAIN;
  // no int NaN exists (quiet_NaN<int>() is 0); an event below every division gets no bin
  int binvalue = -1;
  float value = std::numeric_limits<float>::quiet_NaN();
  for (int i = 0; i < NDIVS; i++)
  {
    if (m_centrality_map[bunch][i] < m_mbd_total_charge)
    {
      binvalue = i + 1;
      value = static_cast<float>(i + 1) / static_cast<float>(NDIVS);
      break;
    }
  }
  if (Verbosity())
  {
    std::cout << " Centile : " << value << std::endl;
    std::cout << " Charge : " << m_mbd_total_charge << std::endl;
  }

  // the (vertex- and run-scaled) charge sum compared with the divisions
  m_central->set_quantity(CentralityInfo::PROP::mbd_NS, m_mbd_total_charge);
  m_central->set_centile(CentralityInfo::PROP::mbd_NS, value);
  if (binvalue >= 0)
  {
    m_central->set_centrality_bin(CentralityInfo::PROP::mbd_NS, binvalue);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << "------------CentralityReco-------------" << std::endl;
  }

  // Get Nodes from the Tree
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  if (!m_mb_info->isAuAuMinimumBias())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (m_bunch_classes)
  {
    // a bunch-class MinimumBiasClassifier sets t0 for every min-bias event; without it (MinimumBiasInfov1, or the
    // classifier run with an older payload) the bunch class is unknown
    if (!std::isfinite(m_mb_info->getT0Corrected()))
    {
      if (!m_warned_no_bunch_class)
      {
        std::cout << PHWHERE << " bunch-class centrality payload, but " << m_mb_info_nodename
                  << " has no bunch class: rerun MinimumBiasClassifier with the same payload. No centrality is filled." << std::endl;
        m_warned_no_bunch_class = true;
      }
      return Fun4AllReturnCodes::EVENT_OK;
    }
    // outside the three bunch windows there is no calibration: centrality stays unset (NaN)
    if (m_mb_info->getBunchClass() == MinimumBiasInfo::BUNCH::NONE)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  // Fill Arrays
  if (FillVars())
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (FillCentralityInfo())
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (Verbosity())
  {
    m_central->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  m_mb_info = findNode::getClass<MinimumBiasInfo>(topNode, m_mb_info_nodename);

  if (!m_mb_info)
  {
    std::cout << "no mb_info node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_central = findNode::getClass<CentralityInfo>(topNode, m_centrality_nodename);

  if (!m_central)
  {
    std::cout << "no centrality node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_mbd_container = findNode::getClass<MbdPmtContainer>(topNode, m_mbd_pmt_nodename);

  if (!m_mbd_container)
  {
    std::cout << "no MBD out node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_mbd_out = findNode::getClass<MbdOut>(topNode, m_mbd_out_nodename);
  if (Verbosity())
  {
    std::cout << "Getting MBD Out" << std::endl;
  }

  if (!m_mbd_out)
  {
    std::cout << "no MBD out node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CentralityReco::CreateNodes(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
  }

  PHNodeIterator dstIter(dstNode);

  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!detNode)
  {
    std::cout << PHWHERE << "Detector Node missing, making one" << std::endl;
    detNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(detNode);
  }

  CentralityInfo *central = new CentralityInfov2();

  PHIODataNode<PHObject> *centralityNode = new PHIODataNode<PHObject>(central, m_centrality_nodename, "PHObject");
  detNode->addNode(centralityNode);

  return;
}

float CentralityReco::getVertexScale(int bunch)
{
  float mbd_vertex = m_mbd_out->get_zvtx();
  for (const auto &v_range_scale : m_vertex_scales[bunch])
  {
    auto v_range = v_range_scale.first;
    if (Verbosity())
    {
      std::cout << "vertexrange : " << v_range.first << "-" << v_range.second << std::endl;
    }

    if (mbd_vertex > v_range.first && mbd_vertex <= v_range.second)
    {
      return v_range_scale.second;
    }
  }
  return 0;
}
