#include "MinimumBiasClassifier.h"

#include "MinimumBiasInfov2.h"

#include <zdcinfo/Zdcinfo.h>

#include <cdbobjects/CDBTTree.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

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
#include <utility>  // for pair

namespace
{
  // class-suffixed payload columns, indexed by MinimumBiasInfo::BUNCH
  const std::array<std::string, 3> bunch_suffix{"", "_m5", "_p5"};
}  // namespace

MinimumBiasClassifier::MinimumBiasClassifier(const std::string &name)
  : SubsysReco(name)
{
}

int MinimumBiasClassifier::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  if (m_species == MinimumBiasInfo::SPECIES::OO)
  {
    m_useZDC         = false;
    m_box_cut        = false;
    m_hit_cut        = 1;
    m_max_charge_cut = 400;
    m_mbd_charge_cut = 0.4;
    m_mbd_time_cut   = 20.;
    m_z_vtx_cut      = 150.;
  }
  else if (m_species == MinimumBiasInfo::SPECIES::PP)
  {
    m_useZDC         = false;
    m_box_cut        = false;
    m_hit_cut        = 1;
    m_max_charge_cut = 300;
    m_mbd_charge_cut = 0.4;
    m_mbd_time_cut   = 20.;
  }

  CDBInterface *m_cdb = CDBInterface::instance();

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

  m_zdc_energy_sum.fill(0);
  m_mbd_charge_sum.fill(0);
  m_mbd_hit.fill(0);

  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::ResetEvent(PHCompositeNode * /*unused*/)
{
  m_zdc_energy_sum.fill(0);
  m_mbd_charge_sum.fill(0);
  m_mbd_hit.fill(0);

  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::FillMinimumBiasInfo()
{
  if (Verbosity())
  {
    std::cout << "Getting Vertex" << std::endl;
  }

  // no t0 yet: events rejected before the classification below have no bunch class
  if (m_bunch_classes)
  {
    m_mb_info->setBunchClass(MinimumBiasInfo::BUNCH::NONE);
  }

  // if (!m_global_vertex_map)
  // {
  //   m_mb_info->setIsAuAuMinimumBias(false);
  //   return Fun4AllReturnCodes::EVENT_OK;
  // }

  // if (m_global_vertex_map->empty())
  // {
  //   m_mb_info->setIsAuAuMinimumBias(false);
  //   return Fun4AllReturnCodes::EVENT_OK;
  // }

  if (m_global_vertex_map->empty())
  {
    m_mb_info->setIsAuAuMinimumBias(false);
    if (m_abortEvents)
    {
      return 1;
    }
    return 0;
  }

  GlobalVertex *vtx = m_global_vertex_map->begin()->second;
  if (m_bunch_classes)
  {
    // the calibration uses the MBD vertex (z and t0)
    std::vector<GlobalVertex *> mbd_vtxs = m_global_vertex_map->get_gvtxs_with_type({GlobalVertex::VTXTYPE::MBD});
    vtx = mbd_vtxs.empty() ? nullptr : mbd_vtxs.front();
  }
  if (!vtx)
  {
    m_mb_info->setIsAuAuMinimumBias(false);
    if (m_abortEvents)
    {
      return 1;
    }
    return 0;
  }

  if (!vtx->isValid())
  {
    m_mb_info->setIsAuAuMinimumBias(false);
    if (m_abortEvents)
    {
      return 1;
    }
    return 0;
  }

  bool minbiascheck = true;

  m_vertex = vtx->get_z();

  int bunch = MinimumBiasInfo::BUNCH::MAIN;
  if (m_bunch_classes)
  {
    const float t0_corrected = vtx->get_t() - m_t0_offset;
    bunch = bunchClass(t0_corrected);
    m_mb_info->setT0Corrected(t0_corrected);
  }
  m_mb_info->setBunchClass(bunch);
  // the bunch windows are the t0 cut: an event outside all three is not min-bias (and has no centrality calibration;
  // its charge sum below uses the main-bunch scales)
  if (bunch == MinimumBiasInfo::BUNCH::NONE)
  {
    minbiascheck = false;
  }
  const int scale_class = (bunch == MinimumBiasInfo::BUNCH::NONE ? MinimumBiasInfo::BUNCH::MAIN : bunch);

  m_vertex_scale = getVertexScale(scale_class);

  if (Verbosity())
  {
    std::cout << "Getting ZDC" << std::endl;
  }
  if (!m_issim && m_useZDC)
  {
    if (!m_zdcinfo)
    {
      m_mb_info->setIsAuAuMinimumBias(false);
      if (m_abortEvents)
      {
        return 1;
      }
      return 0;
    }
  }
  //  Z vertex is within range
  if (std::fabs(m_vertex) > m_z_vtx_cut && minbiascheck)
  {
    minbiascheck = false;
  }
  if (Verbosity())
  {
    std::cout << "Calculating" << std::endl;
  }

  // calculate charge sum and n hit
  for (int i = 0; i < 128; i++)
  {
    m_mbd_pmt = m_mbd_container->get_pmt(i);
    short side = i / 64;
    bool pass = passesHitCut(m_mbd_pmt);
    if (!pass)
    {
      continue;
    }
    m_mbd_hit[side]++;
    m_mbd_charge_sum[side] += m_mbd_pmt->get_q() * m_vertex_scale * m_centrality_scale[scale_class];
  }
  if (Verbosity())
  {
    std::cout << m_mbd_charge_sum[0] << " " << m_mbd_charge_sum[1] << std::endl;
  }

  // MBD Background cut
  if (m_box_cut && !m_bunch_classes && m_mbd_charge_sum[1] < m_mbd_north_cut && m_mbd_charge_sum[0] > m_mbd_south_cut && minbiascheck)
  {
    minbiascheck = false;
    //    m_mb_info->setIsAuAuMinimumBias(false);
    // return Fun4AllReturnCodes::EVENT_OK;
  }

  // Mbd two hit requirement and ZDC energy sum coincidence requirement
  for (int iside = 0; iside < 2; iside++)
  {
    if (m_mbd_hit[iside] < m_hit_cut && minbiascheck)
    {
      minbiascheck = false;
      // m_mb_info->setIsAuAuMinimumBias(false);
      // return Fun4AllReturnCodes::EVENT_OK;
    }
    if (!m_issim && m_useZDC)
    {
      if (m_zdcinfo->get_zdc_energy(iside) <= (m_bunch_classes ? m_bunch_zdc_cut : m_zdc_cut) && minbiascheck)
      {
        minbiascheck = false;
        // m_mb_info->setIsAuAuMinimumBias(false);
        // return Fun4AllReturnCodes::EVENT_OK;
      }
    }
  }
  if ((m_mbd_charge_sum[0] + m_mbd_charge_sum[1]) > m_max_charge_cut && minbiascheck)
  {
    minbiascheck = false;
    // m_mb_info->setIsAuAuMinimumBias(false);
    // return Fun4AllReturnCodes::EVENT_OK;
  }

  if (m_species == MinimumBiasInfo::SPECIES::OO && m_reject_pileup && (m_mbd_charge_sum[0] + m_mbd_charge_sum[1]) > m_pileup_charge_cut && minbiascheck)
  {
    minbiascheck = false;
  }

  m_mb_info->setIsAuAuMinimumBias(minbiascheck);
  if (!minbiascheck && m_abortEvents)
  {
    return 1;
  }
  return 0;
}

int MinimumBiasClassifier::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << "------------MinimumBiasClassifier-------------" << std::endl;
  }

  // Get Nodes from the Tree
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  if (FillMinimumBiasInfo())
  {
    if (Verbosity())
    {
      std::cout << "MinimumBiasClassifier::process_event Aborting Event - not minbias" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  m_mb_info = findNode::getClass<MinimumBiasInfo>(topNode, m_mb_info_nodename);

  if (!m_mb_info)
  {
    std::cout << "no minimum bias node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_mbd_container = findNode::getClass<MbdPmtContainer>(topNode, m_mbd_pmt_nodename);
  if (Verbosity())
  {
    std::cout << "Getting MBD Tubes" << std::endl;
  }

  if (!m_mbd_container)
  {
    std::cout << "no MBD out node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_issim && m_useZDC)
  {
    m_zdcinfo = findNode::getClass<Zdcinfo>(topNode, m_zdc_info_nodename);
    if (Verbosity())
    {
      std::cout << "Getting ZDC Info" << std::endl;
    }

    if (!m_zdcinfo)
    {
      std::cout << "no zdc towers node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  if (Verbosity())
  {
    std::cout << "Getting Vertex Map" << std::endl;
  }

  m_global_vertex_map = findNode::getClass<GlobalVertexMap>(topNode, m_global_vertex_nodename);

  if (!m_global_vertex_map)
  {
    std::cout << "no vertex map node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void MinimumBiasClassifier::CreateNodes(PHCompositeNode *topNode)
{
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
    detNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(detNode);
  }

  std::string nodename = m_mb_info_nodename;
  MinimumBiasInfo *mb = new MinimumBiasInfov2();

  PHIODataNode<PHObject> *mbNode = new PHIODataNode<PHObject>(mb, nodename, "PHObject");
  detNode->addNode(mbNode);

  return;
}

bool MinimumBiasClassifier::passesHitCut(MbdPmtHit *hit) const
{
  if (m_bunch_classes)
  {
    // hit times relative to the run's main-bunch offset, wide enough for all three bunch windows
    return std::fabs(hit->get_time() - m_t0_offset) < m_hit_time_window && hit->get_q() > m_bunch_charge_cut;
  }
  if (fabs(hit->get_time()) >= m_mbd_time_cut)
  {
    return false;
  }
  if (hit->get_q() <= m_mbd_charge_cut)
  {
    return false;
  }

  return true;
}
int MinimumBiasClassifier::Download_centralityScale(const std::string &dbfile)
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
    // bunch-class payload: t0 offset, bunch spacings and per-class scales in the same file (Au+Au data only)
    if (m_use_bunch_classes && !m_issim && m_species == MinimumBiasInfo::SPECIES::AUAU)
    {
      const double t0_offset = cdbttree->GetDoubleValue(0, "t0offset");
      if (std::isfinite(t0_offset))
      {
        m_bunch_classes = true;
        m_t0_offset = t0_offset;
        m_spacing_m = cdbttree->GetDoubleValue(0, "spacing_m");
        m_spacing_p = cdbttree->GetDoubleValue(0, "spacing_p");
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

int MinimumBiasClassifier::Download_centralityVertexScales(const std::string &dbfile)
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
int MinimumBiasClassifier::bunchClass(float t0_corrected) const
{
  // windows centred on each bunch, edges halfway between neighbouring bunches; the inner edges belong to the main
  // bunch (same as the calibration); NaN falls through to NONE
  const double hm = 0.5 * m_spacing_m;
  const double hp = 0.5 * m_spacing_p;
  if (t0_corrected >= -hm && t0_corrected <= hp)
  {
    return MinimumBiasInfo::BUNCH::MAIN;
  }
  if (t0_corrected >= -3 * hm && t0_corrected < -hm)
  {
    return MinimumBiasInfo::BUNCH::M5;
  }
  if (t0_corrected > hp && t0_corrected <= 3 * hp)
  {
    return MinimumBiasInfo::BUNCH::P5;
  }
  return MinimumBiasInfo::BUNCH::NONE;
}

float MinimumBiasClassifier::getVertexScale(int bunch) const
{
  for (const auto &v_range_scale : m_vertex_scales[bunch])
  {
    auto v_range = v_range_scale.first;
    if (m_vertex > v_range.first && m_vertex <= v_range.second)
    {
      return v_range_scale.second;
    }
  }
  return 0;
}
