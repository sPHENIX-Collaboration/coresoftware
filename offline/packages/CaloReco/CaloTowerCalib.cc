#include "CaloTowerCalib.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>
#include <phool/PHCompositeNode.h>

#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfov1.h>

//____________________________________________________________________________..
CaloTowerCalib::CaloTowerCalib(const std::string &name)
  : SubsysReco(name)
  , m_dettype(CaloTowerCalib::HCALOUT)
  , m_detector("HCALOUT")
  , m_DETECTOR(TowerInfoContainer::HCAL)
{
  std::cout << "CaloTowerCalib::CaloTowerCalib(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
CaloTowerCalib::~CaloTowerCalib()
{
  std::cout << "CaloTowerCalib::~CaloTowerCalib() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int CaloTowerCalib::InitRun(PHCompositeNode *topNode)
{
  if (m_dettype == CaloTowerCalib::CEMC)
  {
    m_detector = "CEMC";
    m_DETECTOR = TowerInfoContainer::EMCAL;
    m_fieldname = "cemc_abscalib_mip";
    cdbttree = new CDBTTree("/sphenix/user/shuhangli/DB/cemcDB.root");
  }
  else if (m_dettype == CaloTowerCalib::HCALIN)
  {
    m_detector = "HCALIN";
    m_DETECTOR = TowerInfoContainer::HCAL;
    m_fieldname = "ohcal_abscalib_mip";
    cdbttree = new CDBTTree("/sphenix/user/shuhangli/DB/hcalDB.root");
  }
  else if (m_dettype == CaloTowerCalib::HCALOUT)
  {
    m_detector = "HCALOUT";
    m_DETECTOR = TowerInfoContainer::HCAL;
    m_fieldname = "ohcal_abscalib_mip";
    cdbttree = new CDBTTree("/sphenix/user/shuhangli/DB/hcalDB.root");
  }
  else if (m_dettype == CaloTowerCalib::EPD)
  {
    m_detector = "EPD";
    m_DETECTOR = TowerInfoContainer::SEPD;
  }

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode",
                                                           "DST"));
  if (!dstNode)
  {
    std::cout << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  try
  {
    CreateNodeTree(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  topNode->print();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloTowerCalib::process_event(PHCompositeNode * /*topNode*/)
{
  /*
  int n_channels = _raw_towers->size();

  for (int i = 0; i < n_channels; i++)
  {
    TowerInfov1 *caloinfo_raw = _raw_towers->at(i);

    float raw_amplitude = caloinfo_raw->get_energy();

    TowerInfov1 *caloinfo_calib = new TowerInfov1(*caloinfo_raw);

    int key = _raw_towers->encode_key(i);

    float calibconst = cdbttree->GetFloatValue(key, m_fieldname);

    caloinfo_calib->set_energy(raw_amplitude * calibconst);

    _calib_towers->add(caloinfo_calib, i);
  }
  */
  TowerInfoContainerv1::Range begin_end = _raw_towers->getTowers();
  TowerInfoContainerv1::Iterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
    {
      unsigned int key = rtiter->first;
      TowerInfo *caloinfo_raw = rtiter->second;
      float raw_amplitude = caloinfo_raw->get_energy();

      float calibconst = cdbttree->GetFloatValue(key, m_fieldname);

      unsigned int channel = _calib_towers->decode_key(key);
      
      _calib_towers->at(channel)->set_energy(raw_amplitude * calibconst);
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTowerCalib::CreateNodeTree(PHCompositeNode *topNode)
{
  std::cout << "creating node" << std::endl;
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in RawTowerCalibration::CreateNodes");
  }

  // towers
  std::string RawTowerNodeName = "TOWERS_" + m_detector;
  _raw_towers = findNode::getClass<TowerInfoContainerv1>(dstNode,
                                                         RawTowerNodeName);
  if (!_raw_towers)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << " " << RawTowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + RawTowerNodeName + " node in RawTowerCalibration::CreateNodes");
  }

  std::string CalibTowerNodeName = "TOWERS_Calib_" + m_detector;
  _calib_towers = findNode::getClass<TowerInfoContainerv1>(dstNode,
                                                           CalibTowerNodeName);
  if (!_calib_towers)
  {
    _calib_towers = new TowerInfoContainerv1(m_DETECTOR);

    PHIODataNode<PHObject> *calibtowerNode = new PHIODataNode<PHObject>(_calib_towers, CalibTowerNodeName, "PHObject");
    std::cout << "adding calib tower" << std::endl;
    dstNode->addNode(calibtowerNode);
  }

  return;
}
