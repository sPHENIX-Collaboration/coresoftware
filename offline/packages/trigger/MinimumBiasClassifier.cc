#include "MinimumBiasClassifier.h"

#include "MinimumBiasInfov1.h"

<<<<<<< HEAD
#include <mbd/MbdDefs.h>
#include <mbd/MbdOutV1.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
=======
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <mbd/MbdOut.h>

>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
<<<<<<< HEAD
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <phool/phool.h>
=======
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>
#include <iostream>
#include <map>  // for _Rb_tree_iterator
#include <string>
#include <utility>  // for pair
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

MinimumBiasClassifier::MinimumBiasClassifier(const std::string &name)
  : SubsysReco(name)
{
<<<<<<< HEAD

}

MinimumBiasClassifier::~MinimumBiasClassifier()
{

}

int MinimumBiasClassifier::Init(PHCompositeNode * /*unused*/)
{

  return 0;
=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
}

int MinimumBiasClassifier::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  CreateNodes(topNode);
<<<<<<< HEAD
  return Fun4AllReturnCodes::EVENT_OK;;

}

int MinimumBiasClassifier::ResetEvent(PHCompositeNode *)
=======
  return Fun4AllReturnCodes::EVENT_OK;
  ;
}

int MinimumBiasClassifier::ResetEvent(PHCompositeNode * /*unused*/)
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
<<<<<<< HEAD

  for (int i  =0; i < 2; i++) {
    _zdc_energy_sum[i] = 0.;
    _mbd_charge_sum[i] = 0.;
    _mbd_tubes_hit[i] = 0;
  }

  _z_vertex = 999.99;  
  return Fun4AllReturnCodes::EVENT_OK;;
=======
  _zdc_energy_sum.fill(0);
  _mbd_charge_sum.fill(0);
  _mbd_tubes_hit.fill(0);

  _z_vertex = std::numeric_limits<float>::quiet_NaN();
  _energy = std::numeric_limits<float>::quiet_NaN();
  return Fun4AllReturnCodes::EVENT_OK;
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
}

int MinimumBiasClassifier::FillVars()
{
<<<<<<< HEAD

  if (Verbosity()>1)
=======
  if (Verbosity() > 1)
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  //  if (_global_vertex_map->empty()) return Fun4AllReturnCodes::ABORTEVENT;

<<<<<<< HEAD

=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  GlobalVertex *vtx = _global_vertex_map->begin()->second;

  _z_vertex = vtx->get_z();

  for (int i = 0; i < 2; i++)
<<<<<<< HEAD
    {
      _mbd_charge_sum[i] = _mbd_out->get_q(i);
      _mbd_tubes_hit[i] = _mbd_out->get_q(i);
    }
=======
  {
    _mbd_charge_sum[i] = _mbd_out->get_q(i);
    _mbd_tubes_hit[i] = _mbd_out->get_q(i);
  }
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  if (Verbosity())
  {
    std::cout << "      North: " << _mbd_charge_sum[1] << std::endl;
    std::cout << "      South: " << _mbd_charge_sum[0] << std::endl;
  }

<<<<<<< HEAD
  int j  = 0;
=======
  int j = 0;
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  for (unsigned int i = 0; i < _towers_zdc->size(); i++)
  {
    _tmp_tower = _towers_zdc->get_tower_at_channel(i);
    _energy = _tmp_tower->get_energy();
<<<<<<< HEAD
    if (_energy!=0) {
     
      _zdc_energy_sum[j/3] += _energy;
=======
    if (_energy != 0)
    {
      _zdc_energy_sum[j / 3] += _energy;
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
      j++;
    }
  }
  if (Verbosity())
  {
    std::cout << "      North: " << _zdc_energy_sum[1] << std::endl;
    std::cout << "      South: " << _zdc_energy_sum[0] << std::endl;
  }

<<<<<<< HEAD

  return Fun4AllReturnCodes::EVENT_OK;
}


int MinimumBiasClassifier::FillMinimumBiasInfo()
{
  // Fill is minbias
  

  if (Verbosity()>1)
=======
  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::FillMinimumBiasInfo()
{
  // Fill is minbias

  if (Verbosity() > 1)
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

<<<<<<< HEAD

  bool is_it_min_bias = true;

  // 
  if (fabs(_z_vertex) > _z_vtx_cut) is_it_min_bias = false;


  for (int i = 0; i < 2; i++)
    {
      if (_mbd_tubes_hit[i] < _mbd_tube_cut) is_it_min_bias = false;
      if (_zdc_energy_sum[i] < _zdc_cut) is_it_min_bias = false;
    }

  if (_mbd_charge_sum[1] < _mbd_north_cut && _mbd_charge_sum[0] > _mbd_south_cut) is_it_min_bias = false;
=======
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
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  _mb_info->setIsAuAuMinimumBias(is_it_min_bias);

  return Fun4AllReturnCodes::EVENT_OK;
<<<<<<< HEAD

=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
}

int MinimumBiasClassifier::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << "------------MinimumBiasClassifier-------------" << std::endl;
  }

<<<<<<< HEAD

=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  // Get Nodes from the Tree
  if (GetNodes(topNode))
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Fill Arrays
  if (FillVars())
  {
<<<<<<< HEAD
    std::cout << __LINE__<<std::endl;
=======
    std::cout << __LINE__ << std::endl;
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (FillMinimumBiasInfo())
<<<<<<< HEAD
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  if (Verbosity())
    {
      _mb_info->identify();
    }
=======
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (Verbosity())
  {
    _mb_info->identify();
  }
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  return Fun4AllReturnCodes::EVENT_OK;
}

int MinimumBiasClassifier::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  _mb_info = findNode::getClass<MinimumBiasInfov1>(topNode, "MinimumBiasInfo");
<<<<<<< HEAD
  
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
=======

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
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  _towers_zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_ZDC");

  if (!_towers_zdc)
<<<<<<< HEAD
    {
      std::cout << "no zdc towers node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
=======
  {
    std::cout << "no zdc towers node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac

  _global_vertex_map = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");

  if (!_global_vertex_map)
<<<<<<< HEAD
    {
      std::cout << "no vertex map node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  return Fun4AllReturnCodes::EVENT_OK;

=======
  {
    std::cout << "no vertex map node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
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

<<<<<<< HEAD
  MinimumBiasInfov1 *mb = new MinimumBiasInfov1();
  
=======
  MinimumBiasInfo *mb = new MinimumBiasInfov1();

>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
  PHIODataNode<PHObject> *mbNode = new PHIODataNode<PHObject>(mb, "MinimumBiasInfo", "PHObject");
  detNode->addNode(mbNode);

  return;
}
<<<<<<< HEAD

int MinimumBiasClassifier::End(PHCompositeNode * /* topNode*/)
{

  return 0;
}
=======
>>>>>>> b183955abcd8650f7a6403b3814a7542343be7ac
