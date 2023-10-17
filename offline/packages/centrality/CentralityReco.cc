#include "CentralityReco.h"

#include "CentralityInfov2.h"

#include <bbc/BbcDefs.h>
#include "bbc/BbcOutV1.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

CentralityReco::CentralityReco(const std::string &name)
  : SubsysReco(name)
{

  const int centrality_map[20] = {1999, 1499, 1291, 1102, 937, 790, 660, 547, 449, 363, 289, 227, 174, 130, 94, 66, 45, 0, 0, 0};
  for (int i = 0; i < 20; i++)
  {
    _centrality_map[i] = centrality_map[i];
  }

}

CentralityReco::~CentralityReco()
{

}

int CentralityReco::Init(PHCompositeNode * /*unused*/)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  CreateNodes(topNode);
  return 0;
}

void CentralityReco::ResetVars()
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  _tubes_hit_s = 0;
  _tubes_hit_n = 0;
  _mbd_charge_sum = 0.;
  _mbd_charge_sum_n = 0.;
  _mbd_charge_sum_s = 0.;

  return;
}

int CentralityReco::FillVars()
{

  if (Verbosity()>1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  _mbd_charge_sum_s = _bbc_out->get_q(0);

  _mbd_charge_sum_n = _bbc_out->get_q(1);

  _mbd_charge_sum = _mbd_charge_sum_n + _mbd_charge_sum_s;

  _tubes_hit_s = _bbc_out->get_npmt(0);
  _tubes_hit_n = _bbc_out->get_npmt(1);


  if (Verbosity())
  {
    std::cout << "  MBD sum = " <<_mbd_charge_sum<<std::endl;
    std::cout << "      North: " << _mbd_charge_sum_n << std::endl;
    std::cout << "      South: " << _mbd_charge_sum_s << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int CentralityReco::FillCentralityInfo()
{
  // Fill is minbias
  

  if (Verbosity()>1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }


  float value = -999.99;
  for (int i = 0; i < 20; i++)
  {
    if (_centrality_map[i] < _mbd_charge_sum)
    {
      value = 0.05 * i;
      break;
    }
  }

  _central->set_centile(CentralityInfo::PROP::mbd_NS, value);

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

  // Reset Arrays
  ResetVars();

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
      _central->identify();
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
  }

  _central = findNode::getClass<CentralityInfov2>(topNode, "CentralityInfo");
  
  if (!_central)
    {
      std::cout << "no centrality node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    
  
  _bbc_out = findNode::getClass<BbcOutV1>(topNode, "BbcOut");
  
  if (!_bbc_out)
    {
      std::cout << "no BBC out node " << std::endl;
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

  CentralityInfov2 *central = new CentralityInfov2();
  
  PHIODataNode<PHObject> *centralityNode = new PHIODataNode<PHObject>(central, "CentralityInfo", "PHObject");
  detNode->addNode(centralityNode);

  return;
}

int CentralityReco::End(PHCompositeNode * /* topNode*/)
{

  return 0;
}
