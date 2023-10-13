#include "CentralityReco.h"

#include "CentralityInfov2.h"

#include <bbc/BbcDefs.h>
#include <bbc/BbcVertexMapv1.h>
#include <bbc/BbcVertexv2.h>
#include <bbc/BbcPmtInfoV1.h>
#include <bbc/BbcPmtInfoContainerV1.h>
#include "bbc/BbcOutV1.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

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
  _zdc_gain_factors[0] = 1.37;
  _zdc_gain_factors[1] = 0.64;
  _zdc_gain_factors[2] = 0.44;
  _zdc_gain_factors[3] = 1.39;
  _zdc_gain_factors[4] = 0.78;
  _zdc_gain_factors[5] = 0.29;

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
  _quality = 0;
  _isMinBias = 0;
  _tubes_hit_s = 0;
  _tubes_hit_n = 0;
  _mbd_charge_sum = 0.;
  _mbd_charge_sum_n = 0.;
  _mbd_charge_sum_s = 0.;
  _zdc_energy_sum = 0.;
  _zdc_energy_sum_n = 0.;
  _zdc_energy_sum_s = 0.;

  for (int i = 0; i < 6; i++)
  {
    m_zdc_energy_low[i] = 0;
  }
  return;
}

int CentralityReco::FillVars()
{
  unsigned int size;
  if (Verbosity()>1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  _mbd_charge_sum_s = _bbc_out->get_q(0);
  _mbd_charge_sum_n = _bbc_out->get_q(1);

  _mbd_charge_sum = _mbd_charge_sum_n + _mbd_charge_sum_s;

  _tubes_hit_s = _bbc_out->get_npmt(0);
  _tubes_hit_n = _bbc_out->get_npmt(1);

  _isMinBias = _tubes_hit_s >=2 && _tubes_hit_n >= 2;

  size = _towers_zdc->size();

  if (Verbosity())
  {
    std::cout << "  MBD sum = " <<_mbd_charge_sum<<std::endl;
    std::cout << "      North: " << _mbd_charge_sum_n << std::endl;
    std::cout << "      South: " << _mbd_charge_sum_s << std::endl;
  }
  if (_use_ZDC)
    {
      for (unsigned int i = 0; i < size; i++)
	{

	  _tmp_tower = _towers_zdc->get_tower_at_channel(i);
	  _energy = _tmp_tower->get_energy();

	  if (!(i % 2))
	    {
	      if ((i / 2) % 4 == 3)
		{
		  m_zdc_sum_low[i / 8] = _energy;
		}
	      else
		{
		  m_zdc_energy_low[i / 2] = _zdc_gain_factors[i / 2] * _energy;
		}
	    }
	  else
	    {
	      if ((i / 2) % 4 == 3)
		{
		  m_zdc_sum_high[i / 8] = _energy;
		}
	      else
		{
		  m_zdc_energy_high[i / 2] = _zdc_gain_factors[i / 2] * _energy;
		}
	    }
	}

  
      if (Verbosity() > 5)
	{
	  std::cout << "--------- ZDC data: ----------" << std::endl;
	  std::cout << "South:" << std::endl;
	  for (int i = 0; i < 3; i++)
	    {
	      std::cout << i << " : " << m_zdc_energy_low[i] << " (" << m_zdc_energy_high[i] << ") " << std::endl;
	    }
	  std::cout << "Sum : " << m_zdc_sum_low[0] << " (" << m_zdc_sum_high[0] << ") " << std::endl;
	  std::cout << "North:" << std::endl;
	  for (int i = 0; i < 3; i++)
	    {
	      std::cout << i << " : " << m_zdc_energy_low[i + 3] << " (" << m_zdc_energy_high[i + 3] << ") " << std::endl;
	    }
	  std::cout << "Sum : " << m_zdc_sum_low[1] << " (" << m_zdc_sum_high[1] << ") " << std::endl;
	}
    }

  if (Verbosity() > 5)
    {
    std::cout << "--------- MBD data: ----------" << std::endl;
    std::cout << "South:" << std::endl;
    for (int i = 0; i < 64; i++)
    {
      std::cout << m_mbd_channel[i] << " : " << m_mbd_charge[i] << "  " << m_mbd_time[i] << ") " << std::endl;
    }
    std::cout << "North:" << std::endl;
    for (int i = 64; i < 128; i++)
    {
      std::cout << m_mbd_channel[i] << " : " << m_mbd_charge[i] << "  " << m_mbd_time[i] << ") " << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int CentralityReco::CheckZDC()
{
  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }

  for (unsigned int i = 0; i < 6; i++)
  {
    if (i / 3)
    {
      _zdc_energy_sum_n += m_zdc_energy_low[i];
      _zdc_energy_sum += m_zdc_energy_low[i];
    }
    else
    {
      _zdc_energy_sum_s += m_zdc_energy_low[i];
      _zdc_energy_sum += m_zdc_energy_low[i];
    }
  }
  if (Verbosity())
  {
    std::cout << "  ZDC sums: " <<std::endl;
    std::cout << "      North: " << _zdc_energy_sum_n << std::endl;
    std::cout << "      South: " << _zdc_energy_sum_s << std::endl;
  }

  _zdc_check = (_zdc_energy_sum_n > _zdc_energy_threshold) && (_zdc_energy_sum_s > _zdc_energy_threshold);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CentralityReco::FillCentralityInfo()
{
  // Fill is minbias
  
  bool final_check = _isMinBias;
  if (_use_ZDC) final_check |= _zdc_check;
  _central->setMinBias(final_check);

  float value = 0;
  for (int i = 0; i < 20; i++)
  {
    if (_centrality_map[i] < _mbd_charge_sum)
    {
      value = 0.05 * i;
      break;
    }
  }
  if (!value)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
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

  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << " :: " << __LINE__ << std::endl;
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

  // Calculate the charge
  if (GetMBDCharge())
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (_use_ZDC)
    {
      // Check the ZDC coincidence
      if (CheckZDC())
	{
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
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

  if (_use_ZDC)
    {
      _towers_zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");
      
      if (!_towers_zdc)
	{
	  std::cout << "no zdc towers node " << std::endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
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
