// $Id: $

/*!
 * \file DeadHotMapLoader.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "DeadHotMapLoader.h"

#include <calobase/RawTowerDeadMap.h>
#include <calobase/RawTowerDeadMapv1.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/TowerInfoDefs.h>

#include <cdbobjects/CDBTTree.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>                 // for PHWHERE

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

DeadHotMapLoader::DeadHotMapLoader(const std::string &detector)
  : SubsysReco("DeadHotMapLoader_" + detector)
  , m_detector(detector)
{
}

int DeadHotMapLoader::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << Name() << "::" << m_detector << "::" << __PRETTY_FUNCTION__
              << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find Run node in RawTowerCalibration::CreateNodes");
  }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(runNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_detector);
    runNode->addNode(DetNode);
  }

  // Be careful as a previous calibrator may have been registered for this detector
  std::string deadMapName = "DEADMAP_" + m_detector;
  RawTowerDeadMap *m_deadmap = findNode::getClass<RawTowerDeadMapv1>(DetNode, deadMapName);
  if (!m_deadmap)
  {
    const RawTowerDefs::CalorimeterId caloid = RawTowerDefs::convert_name_to_caloid(m_detector);

    m_deadmap = new RawTowerDeadMapv1(caloid);
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_deadmap, deadMapName, "PHObject");
    DetNode->addNode(towerNode);
  }

  // assert(m_deadmap);

  std::string url = CDBInterface::instance()->getUrl(m_detector + "_BadTowerMap");
  if(url.empty())
  {
    std::cout << PHWHERE << " Could not get Dead Map for CDB. Detector: " << m_detector << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_CDBTTree = new CDBTTree(url);

  if(!m_CDBTTree)
  {
    std::cout << "No CDB TTree found from url " << url << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  int etabins;
  int phibins;

  if (m_detector.c_str()[0] == 'H')
    {//HCal towers
      etabins = 24;
      phibins = 64;

     
      for(int i = 0; i < etabins*phibins; i++)
	{
	  int isDead = m_CDBTTree->GetIntValue(i,"status");
	  if(isDead > 0)
	    {
	      unsigned int key = TowerInfoDefs::encode_hcal(i);
	      int ieta = TowerInfoDefs::getCaloTowerEtaBin(key);
	      int iphi = TowerInfoDefs::getCaloTowerPhiBin(key);
	      m_deadmap->addDeadTower(ieta, iphi);
	    }
	}
    }
  
  else if(m_detector.c_str()[0] == 'C')
    {
      etabins = 96;
      phibins = 256;

      for(int i = 0; i < 96*256; i++)
	{
	  int isDead = m_CDBTTree->GetIntValue(i,"status");
	  if(isDead > 0)
	    {
	      unsigned int key = TowerInfoDefs::encode_emcal(i);
	      int ieta = TowerInfoDefs::getCaloTowerEtaBin(key);
	      int iphi = TowerInfoDefs::getCaloTowerPhiBin(key);
	      m_deadmap->addDeadTower(ieta, iphi);
	    }
	}
    }
  

  if (Verbosity())
  {
    std::cout << "DeadHotMapLoader::" << m_detector << "::InitRun - loading dead map completed : ";
    m_deadmap->identify();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
