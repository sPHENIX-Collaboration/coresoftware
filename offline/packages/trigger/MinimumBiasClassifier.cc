#include "MinimumBiasClassifier.h"

#include "MinimumBiasInfov1.h"

#include <zdcinfo/Zdcinfo.h>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

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
#include <map>  // for _Rb_tree_iterator
#include <string>
#include <utility>  // for pair

MinimumBiasClassifier::MinimumBiasClassifier(const std::string &name)
  : SubsysReco(name)
  , m_MinBiasParams(name)
{
}
int MinimumBiasClassifier::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  CDBInterface *m_cdb = CDBInterface::instance();

  std::string centscale_url = m_cdb->getUrl("CentralityScale");
  if (m_overwrite_scale)
  {
    centscale_url = m_overwrite_url_scale;
    std::cout << " Overwriting Scale to " << m_overwrite_url_scale << std::endl;
  }

  if (Download_centralityScale(centscale_url))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  std::string vertexscale_url = m_cdb->getUrl("CentralityVertexScale");
  if (m_overwrite_vtx)
  {
    vertexscale_url = m_overwrite_url_vtx;
    std::cout << " Overwriting Vtx to " << m_overwrite_url_vtx << std::endl;
  }

  if (Download_centralityVertexScales(vertexscale_url))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CreateNodes(topNode);

  // Create Space on NodeTree to save Minimum Bias Params
  PHNodeIterator parIter(topNode);
  m_parNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", "PAR"));
  if (!m_parNode)
  {
    std::cout << "No RUN node found; cannot create PHParameters. Aborting run!";
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_MinBiasParams.SaveToNodeTree(m_parNode, "MinBiasParams");

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

  // set defaults
  m_MinBiasParams.set_int_param("minbias_background_cut_fail", false);
  m_MinBiasParams.set_int_param("minbias_two_hit_min_fail", false);
  m_MinBiasParams.set_int_param("minbias_zdc_energy_min_fail", false);
  m_MinBiasParams.set_int_param("minbias_mbd_total_energy_max_fail", false);
  m_MinBiasParams.UpdateNodeTree(m_parNode, "MinBiasParams");

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
    return Fun4AllReturnCodes::EVENT_OK;
  }

  GlobalVertex *vtx = m_global_vertex_map->begin()->second;
  if (!vtx)
  {
    m_mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (!vtx->isValid())
  {
    m_mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  bool minbiascheck = true;
  ;

  m_vertex = vtx->get_z();

  m_vertex_scale = getVertexScale();

  if (Verbosity())
  {
    std::cout << "Getting ZDC" << std::endl;
  }
  if (!m_issim)
  {
    if (!m_zdcinfo)
    {
      m_mb_info->setIsAuAuMinimumBias(false);
      return Fun4AllReturnCodes::EVENT_OK;
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
    m_mbd_charge_sum[side] += m_mbd_pmt->get_q() * m_vertex_scale * m_centrality_scale;
  }

  m_MinBiasParams.set_double_param("minbias_mbd_total_charge_south", m_mbd_charge_sum[0]);
  m_MinBiasParams.set_double_param("minbias_mbd_total_charge_north", m_mbd_charge_sum[1]);
  m_MinBiasParams.set_double_param("minbias_vertex_scale", m_vertex_scale);
  m_MinBiasParams.set_double_param("minbias_centrality_scale", m_centrality_scale);

  if (Verbosity())
  {
    std::cout << m_mbd_charge_sum[0] << " " << m_mbd_charge_sum[1] << std::endl;
  }

  // MBD Background cut
  if (m_mbd_charge_sum[1] < m_mbd_north_cut && m_mbd_charge_sum[0] > m_mbd_south_cut && minbiascheck)
  {
    minbiascheck = false;
    m_MinBiasParams.set_int_param("minbias_background_cut_fail", true);
    //    m_mb_info->setIsAuAuMinimumBias(false);
    // return Fun4AllReturnCodes::EVENT_OK;
  }

  // Mbd two hit requirement and ZDC energy sum coincidence requirement
  for (int iside = 0; iside < 2; iside++)
  {
    if (m_mbd_hit[iside] < 2 && minbiascheck)
    {
      minbiascheck = false;
      m_MinBiasParams.set_int_param("minbias_two_hit_min_fail", true);
      // m_mb_info->setIsAuAuMinimumBias(false);
      // return Fun4AllReturnCodes::EVENT_OK;
    }
    if (!m_issim)
    {
      if (m_zdcinfo->get_zdc_energy(iside) <= m_zdc_cut && minbiascheck)
      {
        minbiascheck = false;
        m_MinBiasParams.set_int_param("minbias_zdc_energy_min_fail", true);
        // m_mb_info->setIsAuAuMinimumBias(false);
        // return Fun4AllReturnCodes::EVENT_OK;
      }
    }
  }
  if ((m_mbd_charge_sum[0] + m_mbd_charge_sum[1]) > m_mbd_total_charge_cut && minbiascheck)
  {
    minbiascheck = false;
    m_MinBiasParams.set_int_param("minbias_mbd_total_energy_max_fail", true);
    // m_mb_info->setIsAuAuMinimumBias(false);
    // return Fun4AllReturnCodes::EVENT_OK;
  }

  m_MinBiasParams.UpdateNodeTree(m_parNode, "MinBiasParams");
  m_mb_info->setIsAuAuMinimumBias(minbiascheck);

  return Fun4AllReturnCodes::EVENT_OK;
}
int MinimumBiasClassifier::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << "------------MinimumBiasClassifier-------------" << std::endl;
  }

  // Get Nodes from the Tree
  if (GetNodes(topNode))
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (FillMinimumBiasInfo())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  m_mb_info = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");

  if (!m_mb_info)
  {
    std::cout << "no minimum bias node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_mbd_container = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if (Verbosity())
  {
    std::cout << "Getting MBD Tubes" << std::endl;
  }

  if (!m_mbd_container)
  {
    std::cout << "no MBD out node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_issim)
  {
    m_zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
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

  m_global_vertex_map = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

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
    std::cout << PHWHERE << "Detector Node missing, making one" << std::endl;
    detNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(detNode);
  }

  MinimumBiasInfo *mb = new MinimumBiasInfov1();

  PHIODataNode<PHObject> *mbNode = new PHIODataNode<PHObject>(mb, "MinimumBiasInfo", "PHObject");
  detNode->addNode(mbNode);

  return;
}

bool MinimumBiasClassifier::passesHitCut(MbdPmtHit *hit) const
{
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
  m_centrality_scale = 1.00;

  std::filesystem::path dbase_file = dbfile;
  if (dbase_file.extension() == ".root")
  {
    CDBTTree *cdbttree = new CDBTTree(dbase_file);
    cdbttree->LoadCalibrations();
    m_centrality_scale = cdbttree->GetDoubleValue(0, "centralityscale");
    if (Verbosity())
    {
      std::cout << "centscale = " << m_centrality_scale << std::endl;
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

    for (int iv = 0; iv < nvertexbins; iv++)
    {
      float scale = cdbttree->GetDoubleValue(iv, "scale");
      float lowvertex = cdbttree->GetDoubleValue(iv, "low_vertex");
      float highvertex = cdbttree->GetDoubleValue(iv, "high_vertex");
      m_vertex_scales.emplace_back(std::make_pair(lowvertex, highvertex), scale);
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
float MinimumBiasClassifier::getVertexScale()
{
  for (auto v_range_scale : m_vertex_scales)
  {
    auto v_range = v_range_scale.first;
    if (m_vertex > v_range.first && m_vertex <= v_range.second)
    {
      return v_range_scale.second;
    }
  }
  return 0;
}
