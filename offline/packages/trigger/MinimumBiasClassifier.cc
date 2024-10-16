#include "MinimumBiasClassifier.h"

#include "MinimumBiasInfov1.h"

#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <zdcinfo/Zdcinfo.h>

#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>

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
  ;
}

int MinimumBiasClassifier::ResetEvent(PHCompositeNode * /*unused*/)
{
  _zdc_energy_sum.fill(0);
  _mbd_charge_sum.fill(0);

  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::FillMinimumBiasInfo()
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  bool is_it_min_bias = true;

  if (!_global_vertex_map)
  {
    _mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (_global_vertex_map->empty())
  {
    _mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  GlobalVertex *vtx = _global_vertex_map->begin()->second;

  if (!vtx)
  {
    std::cout << "nothing in vertex " << std::endl;

    _mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (!vtx->isValid())
  {
    std::cout << "invalid vertex " << std::endl;

    _mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (!_zdcinfo)
  {
    std::cout << "nothing in zdc " << std::endl;

    _mb_info->setIsAuAuMinimumBias(false);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (std::fabs(vtx->get_z()) > _z_vtx_cut)
  {
    is_it_min_bias = false;
  }

  int npmt = 128;//_pmts_mbd->get_npmt();
  for ( int ipmt = 0; ipmt < npmt; ipmt++)
    {
      _tmp_pmt = _mbd_pmt->get_pmt(ipmt);

      float charge = _tmp_pmt->get_q();
      float time = _tmp_pmt->get_time();
      int side = ipmt/64;
      if (time < time_threshold && charge > charge_threshold)
	{
	  _mbd_charge_sum[side] += charge;
	}
    }

  for (int i = 0; i < 2; i++)
  {
    if (_mbd_charge_sum[i] < _mbd_tube_cut)
    {
      is_it_min_bias = false;
    }
    if (_zdcinfo->get_zdc_energy(i) < _zdc_cut)
    {
      is_it_min_bias = false;
    }
  }

  if (_mbd_charge_sum[1] < _mbd_north_cut && _mbd_charge_sum[0] > _mbd_south_cut)
  {
    is_it_min_bias = false;
  }
  _mb_info->setIsAuAuMinimumBias(is_it_min_bias);
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
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (FillMinimumBiasInfo())
  {
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

  _mb_info = findNode::getClass<MinimumBiasInfo>(topNode, "MinimumBiasInfo");

  if (!_mb_info)
  {
    std::cout << "no minimum bias node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _mbd_pmt = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");

  if (!_mbd_pmt)
  {
    std::cout << "no MBD pmt node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  _zdcinfo = findNode::getClass<Zdcinfo>(topNode, "Zdcinfo");
  if (!_zdcinfo)
  {
    std::cout << "no zdc info node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _global_vertex_map = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

  if (!_global_vertex_map)
  {
    std::cout << "no vertex map node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void MinimumBiasClassifier::CreateNodes(PHCompositeNode *topNode)
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

  MinimumBiasInfo *mb = new MinimumBiasInfov1();

  PHIODataNode<PHObject> *mbNode = new PHIODataNode<PHObject>(mb, "MinimumBiasInfo", "PHObject");
  detNode->addNode(mbNode);

  return;
}
