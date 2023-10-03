#include "CentralityReco.h"

#include "CentralityInfov2.h"
#include "bbc/BbcDefs.h"
#include "bbc/BbcVertexMapv1.h"
#include "bbc/BbcVertexv2.h"

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

  _offset = 28.52;

}

CentralityReco::~CentralityReco()
{

}

int CentralityReco::Init(PHCompositeNode * /*unused*/)
{
  // Histograms

  std::fill_n(gaincorr, 128, 1);
  std::fill_n(tq_t0_offsets, 128, 1);
  const char *gaincalib = getenv("CENTRALITY_GAINCALIB");
  std::string gainfilename;
  if (gaincalib == nullptr)
  {
    const char *offline_main = getenv("OFFLINE_MAIN");
    assert(offline_main);  // make cppcheck happy
    gainfilename = offline_main;
    gainfilename += "/share/centrality/gainfile.calib";
  }
  else
  {
    gainfilename = gaincalib;
  }
  if (!std::filesystem::exists(gainfilename))
  {
    std::cout << PHWHERE << gainfilename << " does not exist" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  std::ifstream gainfile(gainfilename);

  int ch;
  float integ;
  float integerr;
  float peak;
  float peakerr;
  float mbd_width;
  float widtherr;
  float chi2ndf;

  while (gainfile >> ch >> integ >> peak >> mbd_width >> integerr >> peakerr >> widtherr >> chi2ndf)
  {
    gaincorr[ch] = 1.0 / peak;
  }

  gainfile.close();

  const char *bbc_tq_t0calib = getenv("CENTRALITY_BBC_TQ_T0CALIB");
  std::string bbc_tq_t0_filename;
  if (bbc_tq_t0calib == nullptr)
  {
    const char *offline_main = getenv("OFFLINE_MAIN");
    assert(offline_main);  // make cppcheck happy
    bbc_tq_t0_filename = offline_main;
    bbc_tq_t0_filename += "/share/centrality/bbc_tq_t0.calib";
  }
  else
  {
    bbc_tq_t0_filename = bbc_tq_t0calib;
  }
  if (!std::filesystem::exists(bbc_tq_t0_filename))
  {
    std::cout << PHWHERE << bbc_tq_t0_filename << " does not exist" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  std::ifstream tfile(bbc_tq_t0_filename);

  int pmtnum;
  float meanerr;
  float sigma;
  float sigmaerr;
  for (float &tq_t0_offset : tq_t0_offsets)
  {
    tfile >> pmtnum >> tq_t0_offset >> meanerr >> sigma >> sigmaerr;
  }

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
  _tubes_hit[0] = 0;
  _tubes_hit[1] = 0;
  _mbd_charge_sum = 0.;
  _mbd_charge_sum_n = 0.;
  _mbd_charge_sum_s = 0.;
  _zdc_energy_sum = 0.;
  _zdc_energy_sum_n = 0.;
  _zdc_energy_sum_s = 0.;

  for (int i = 0; i < 128; i++)
  {
    m_mbd_charge[i] = 0.;
    m_mbd_time[i] = 0.;
    m_mbd_charge_raw[i] = 0.;
    m_mbd_time_raw[i] = 0.;
    m_mbd_side[i] = 0;
    m_mbd_channel[i] = 0;
  }
  for (int i = 0; i < 6; i++)
  {
    m_zdc_energy_low[i] = 0;
    m_zdc_energy_high[i] = 0;
  }
  for (int i = 0; i < 2; i++)
  {
    m_zdc_sum_low[i] = 0;
    m_zdc_sum_high[i] = 0;
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
  size = _towers_mbd->size();
  
  if (size != 256)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  for (unsigned int i = 0; i < size; i++)
    {
      _tmp_tower = _towers_mbd->get_tower_at_channel(i);
      _key = TowerInfoDefs::encode_mbd(i);
      _type = TowerInfoDefs::get_mbd_type(_key);
      _side = TowerInfoDefs::get_mbd_side(_key);
      _channel = TowerInfoDefs::get_mbd_channel(_key);
      _energy = _tmp_tower->get_energy();
      if (_type)
	{
	  m_mbd_charge_raw[_channel + _side * 64] = _energy;
	  m_mbd_charge[_channel + _side * 64] = gaincorr[_channel + _side * 64] * _energy;
	}
      else
	{
	  m_mbd_time_raw[_channel + _side * 64] = _energy;
	  m_mbd_time[_channel + _side * 64] = (25. - _energy * (9.0 / 5000.) - tq_t0_offsets[_channel + _side * 64]);
	  m_mbd_side[_channel + _side * 64] = _side;
	  m_mbd_channel[_channel + _side * 64] = _channel;
	}
    }

  size = _towers_zdc->size();

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

int CentralityReco::GetMBDCharge()
{

  if (Verbosity() > 1)
  {
    std::cout << __FILE__ << " :: " << __FUNCTION__ << std::endl;
  }
  _tdc[0] = 0.;
  _tdc[1] = 0.;
  unsigned int size = _towers_mbd->size();
  int side;

  for (unsigned int i = 0; i < size / 2; i++)
  {
    side = m_mbd_side[i];
    if (m_mbd_charge[i] <= _mbd_charge_threshold)
    {
      continue;
    }

    _tubes_hit[side]++;
  }

  if ((_tubes_hit[0] >= 2) && (_tubes_hit[1] >= 2))
  {
    _isMinBias = 1;
  }
  else
  {
    _isMinBias = 0;
  }

  for (unsigned int i = 0; i < size / 2; i++)
  {
    side = m_mbd_side[i];

    if (side)
    {
      _mbd_charge_sum_n += m_mbd_charge[i];
    }
    else
    {
      _mbd_charge_sum_s += m_mbd_charge[i];
    }
  }

  _mbd_charge_sum = _mbd_charge_sum_s + _mbd_charge_sum_n;
  
  if (Verbosity())
  {
    std::cout << "  Mbd Charge Sum = " << _mbd_charge_sum << std::endl;
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
  bool final_check = _zdc_check && _isMinBias;
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

  // Check the ZDC coincidence
  if (CheckZDC())
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
    
  
  _towers_mbd = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_MBD");
  
  if (!_towers_mbd)
    {
      std::cout << "no mbd towers node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  _towers_zdc = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_ZDC");

  if (!_towers_zdc)
  {
    std::cout << "no zdc towers node " << std::endl;
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
