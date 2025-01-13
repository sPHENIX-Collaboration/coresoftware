#include "CentralityReco.h"

#include "CentralityInfov2.h"

#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <calotrigger/MinimumBiasInfo.h>

#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>

#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBTTree.h>

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

CentralityReco::CentralityReco(const std::string &name)
  : SubsysReco(name)
{

}

int CentralityReco::InitRun(PHCompositeNode *topNode)
{
  CDBInterface *m_cdb = CDBInterface::instance();

  std::string centdiv_url = m_cdb->getUrl("Centrality");
  if (m_overwrite_divs)
    {
      centdiv_url = m_overwrite_url_divs;
      std::cout << " Overwriting Divs to " << m_overwrite_url_divs << std::endl;
    }

  if (Download_centralityDivisions(centdiv_url))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

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
  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::Download_centralityScale(const std::string &dbfile)
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

int CentralityReco::Download_centralityDivisions(const std::string &dbfile)
{
  m_centrality_map.fill(0);

  std::filesystem::path dbase_file = dbfile;

  if (dbase_file.extension() == ".root")
  {
    CDBTTree *cdbttree = new CDBTTree(dbase_file);
    cdbttree->LoadCalibrations();
    if (Verbosity())
      {
	cdbttree->Print();
      }
    for (int idiv = 0; idiv < NDIVS; idiv++)
    {
      m_centrality_map[idiv] = cdbttree->GetFloatValue(idiv, "centralitydiv");
      if (Verbosity())
      {
        std::cout << "centdiv " << idiv << " : " << m_centrality_map[idiv] << std::endl;
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


  float scale_factor = getVertexScale();

  if (Verbosity())
    {
      std::cout << scale_factor << "*" << m_centrality_scale << std::endl;
    }

  for (int i = 0; i < 128; i++)
    {
      
      m_mbd_hit = m_mbd_container->get_pmt(i);

      if ((m_mbd_hit->get_q()) < mbd_charge_cut)
	{
	  continue;
	}
      if (fabs(m_mbd_hit->get_time()) > mbd_time_cut) 
	{
	  continue;
	}
      m_mbd_total_charge += m_mbd_hit->get_q()*scale_factor*m_centrality_scale;
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

  int binvalue = std::numeric_limits<int>::quiet_NaN();
  float value = std::numeric_limits<float>::quiet_NaN();
  for (int i = 0; i < NDIVS; i++)
  {
    if (m_centrality_map[i] < m_mbd_total_charge)
    {
      binvalue = i + 1;
      value =  static_cast<float>(i + 1)/static_cast<float>(NDIVS);
      break;
    }
  }
  if (Verbosity()) 
    {
      std::cout << " Centile : " << value << std::endl;      
      std::cout << " Charge : " << m_mbd_total_charge << std::endl;
    }

  m_central->set_centile(CentralityInfo::PROP::mbd_NS, value);
  m_central->set_centrality_bin(CentralityInfo::PROP::mbd_NS, binvalue);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << "------------CentralityReco-------------" << std::endl;
  }

  // Get Nodes from the Tree
  if (GetNodes(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!m_mb_info->isAuAuMinimumBias())
    {
      return Fun4AllReturnCodes::EVENT_OK;
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

  m_mb_info = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");

  if (!m_mb_info)
  {
    std::cout << "no mb_info node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_global_vertex_map = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  
  if (!m_global_vertex_map)
    {
    std::cout << "no vertex map node " << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_central = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");

  if (!m_central)
  {
    std::cout << "no centrality node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_mbd_container = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");

  if (!m_mbd_container)
  {
    std::cout << "no MBD out node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }


  m_mbd_out = findNode::getClass<MbdOut>(topNode, "MbdOut");
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

  PHIODataNode<PHObject> *centralityNode = new PHIODataNode<PHObject>(central, "CentralityInfo", "PHObject");
  detNode->addNode(centralityNode);

  return;
}

float CentralityReco::getVertexScale()
{


  float mbd_vertex = m_mbd_out->get_zvtx();

  for (auto v_range_scale : m_vertex_scales)
    {
      auto v_range = v_range_scale.first;
      if (Verbosity())
	{
	  std::cout << "vertexrange : "<<v_range.first<<"-"<<v_range.second << std::endl;
	}

      if (mbd_vertex > v_range.first && mbd_vertex <= v_range.second)
	{
	  return v_range_scale.second;
	}
    }
  return 0;
}

