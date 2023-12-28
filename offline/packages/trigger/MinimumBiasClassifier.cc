#include "MinimumBiasClassifier.h"

#include "MinimumBiasInfov1.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdOut.h>

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
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  _zdc_energy_sum.fill(0);
  _mbd_charge_sum.fill(0);
  _mbd_tubes_hit.fill(0);

  _z_vertex = std::numeric_limits<float>::quiet_NaN();
  _energy = std::numeric_limits<float>::quiet_NaN();
  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::FillVars()
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  //  if (_global_vertex_map->empty()) return Fun4AllReturnCodes::ABORTEVENT;

  GlobalVertex *vtx = _global_vertex_map->begin()->second;

  _z_vertex = vtx->get_z();

  for (int i = 0; i < 2; i++)
  {
    _mbd_charge_sum[i] = _mbd_out->get_q(i);
    _mbd_tubes_hit[i] = _mbd_out->get_q(i);
  }

  if (Verbosity())
  {
    std::cout << "      North: " << _mbd_charge_sum[1] << std::endl;
    std::cout << "      South: " << _mbd_charge_sum[0] << std::endl;
  }

  int j = 0;
  for (unsigned int i = 0; i < _towers_zdc->size(); i++)
  {
    _tmp_tower = _towers_zdc->get_tower_at_channel(i);
    _energy = _tmp_tower->get_energy();
    if (_energy != 0)
    {
      _zdc_energy_sum[j / 3] += _energy;
      j++;
    }
  }
  if (Verbosity())
  {
    std::cout << "      North: " << _zdc_energy_sum[1] << std::endl;
    std::cout << "      South: " << _zdc_energy_sum[0] << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::FillMinimumBiasInfo()
{
  // Fill is minbias

  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  bool is_it_min_bias = true;

  //
  if (std::fabs(_z_vertex) > _z_vtx_cut)
  {
    is_it_min_bias = false;
  }

  for (int i = 0; i < 2; i++)
  {
    if (_mbd_tubes_hit[i] < _mbd_tube_cut)
    {
      is_it_min_bias = false;
    }
    if (_zdc_energy_sum[i] < _zdc_cut)
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

  // Fill Arrays
  if (FillVars())
  {
    std::cout << __LINE__ << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (FillMinimumBiasInfo())
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (Verbosity())
  {
    _mb_info->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  _mb_info = findNode::getClass<MinimumBiasInfov1>(topNode, "MinimumBiasInfo");

  if (!_mb_info)
  {
    std::cout << "no minimum bias node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _mbd_out = findNode::getClass<MbdOut>(topNode, "MbdOut");

  if (!_mbd_out)
  {
    std::cout << "no MBD out node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _towers_zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_ZDC");

  if (!_towers_zdc)
  {
    std::cout << "no zdc towers node " << std::endl;
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
