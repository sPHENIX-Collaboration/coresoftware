#include "MinimumBiasClassifier.h"

#include "MinimumBiasInfov1.h"

#include <zdcinfo/Zdcinfo.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdOut.h>
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
#include <iostream>
#include <map>  // for _Rb_tree_iterator
#include <string>
#include <utility>  // for pair

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

  CreateNodes(topNode);

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

  if (!m_global_vertex_map)
  {
    m_mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;
  }

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
  if (Verbosity())
    {
      std::cout << "Getting ZDC" << std::endl;
    }

  if (!m_zdcinfo)
  {
    m_mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  //  Z vertex is within range
  if (std::fabs(vtx->get_z()) > m_z_vtx_cut)
  {
    m_mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;  
  }
  if (Verbosity())
    {
      std::cout << "Calculating" << std::endl;
    }

  // calculate charge sum and n hit
  for (int i = 0; i < 128; i++)
  {
    m_mbd_pmt = m_mbd_container->get_pmt(i);
    short side = i/64;
    bool pass = passesHitCut(m_mbd_pmt);
    if (!pass)
    {
      continue;
    }
    m_mbd_hit[side]++;
    m_mbd_charge_sum[side] += m_mbd_pmt->get_q();
  }
  if (Verbosity())
    {
      std::cout <<m_mbd_charge_sum[0] << " " << m_mbd_charge_sum[1]<< std::endl;
    }

  // MBD Background cut
  if (m_mbd_charge_sum[1] < m_mbd_north_cut && m_mbd_charge_sum[0] > m_mbd_south_cut)
  {
    m_mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;  
  }

  // Mbd two hit requirement and ZDC energy sum coincidence requirement
  for (int iside = 0; iside < 2; iside++)
    {
      if (m_mbd_hit[iside] < 2)
	{
	  m_mb_info->setIsAuAuMinimumBias(false);
	  return Fun4AllReturnCodes::EVENT_OK;
	}
      if (m_zdcinfo->get_zdc_energy(iside) <= m_zdc_cut)
	{
	  m_mb_info->setIsAuAuMinimumBias(false);
	  return Fun4AllReturnCodes::EVENT_OK;	 
	}
    }

  m_mb_info->setIsAuAuMinimumBias(true);
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
 bool MinimumBiasClassifier::passesHitCut(MbdPmtHit *hit)
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
