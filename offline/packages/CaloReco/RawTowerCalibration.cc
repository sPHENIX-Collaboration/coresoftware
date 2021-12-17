#include "RawTowerCalibration.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv1.h>
#include <calobase/RawTowerv2.h>

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <fstream>
using namespace std;

RawTowerCalibration::RawTowerCalibration(const std::string &name)
  : SubsysReco(name)
  ,  //
  calibfile("empty.txt")
  , _tower_calib_params(name)
{
}

int RawTowerCalibration::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode",
                                                           "DST"));
  if (!dstNode)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerCalibration::process_event(PHCompositeNode */*topNode*/)
{
   if (Verbosity())
   {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << "Process event entered" << std::endl;
   } 

  RawTowerContainer::ConstRange begin_end = m_RawTowerContainer->getTowers();
  RawTowerContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
  {
    const RawTowerDefs::keytype key = rtiter->first;
    
    const RawTower *raw_tower = rtiter->second;
    assert(raw_tower);

    RawTowerGeom *raw_tower_geom = m_RawTowerGeomContainer->get_tower_geometry(raw_tower->get_id());
    assert(raw_tower_geom);

    if (m_TowerType >= 0)
    {
      // Skip towers that don't match the type we are supposed to calibrate
      if (m_TowerType != raw_tower_geom->get_tower_type())
      {
        continue;
      }
    }

    if (m_CalibAlgorithmEnum == kNo_calibration)
    {
      m_CalibTowerContainer->AddTower(key, new RawTowerv2(*raw_tower));
    }
    else if (m_CalibAlgorithmEnum == kSimple_linear_calibration)
    {
      const double raw_energy = raw_tower->get_energy();
      const double calib_energy = (raw_energy - m_PedestalADC) * m_CalibConst_GeV_per_ADC;
      
       RawTower *calib_tower = new RawTowerv2(*raw_tower);
      calib_tower->set_energy(calib_energy);
      m_CalibTowerContainer->AddTower(key, calib_tower);
    }
    else if (m_CalibAlgorithmEnum == kTower_by_tower_calibration)
    {
      RawTowerDefs::CalorimeterId caloid = RawTowerDefs::decode_caloid(key);
      const int eta = raw_tower->get_bineta();
      const int phi = raw_tower->get_binphi();
 
      double tower_by_tower_calib = 1.;
      double recalibrated_e = 0.0;
      std::string re_cal_flag = "empty";
 
     if (caloid == RawTowerDefs::LFHCAL) 
      {
        const int l   = raw_tower->get_binl();
        const string calib_const_name("calib_const_eta" + to_string(eta) + "_phi" + to_string(phi) + "_l" + to_string(l));

        tower_by_tower_calib = _tower_calib_params.get_double_param(calib_const_name);

        if (m_PedestalFromFileFlag == true)
        {
          const string pedstal_name("PedCentral_ADC_eta" + to_string(eta) + "_phi" + to_string(phi)+ "_l" + to_string(l));
          m_PedestalADC = _tower_calib_params.get_double_param(pedstal_name);
        }

        if (m_GeV_per_ADC_FromFileFlag == true)
        {
          const string GeVperADCname("GeVperADC_eta" + to_string(eta) + "_phi" + to_string(phi)+ "_l" + to_string(l));
          m_CalibConst_GeV_per_ADC = _tower_calib_params.get_double_param(GeVperADCname);
        }
      } 
      else
      {     

       double recal_e[24][64] = {{0.0}};
       ifstream calibrate_tower;
       calibrate_tower.open(calibfile);
       if (calibrate_tower.is_open())
        {
            int rows = 0;
            while (!calibrate_tower.eof())
            {
              int etabin = -1;
              int phibin = -1;
              double recal = 0.0;
              calibrate_tower >> etabin >> phibin >> recal;
              recal_e[etabin][phibin] = recal;
              rows++;
              if (rows > 1)
              {
                re_cal_flag = "CALIBMODE";
              }
            }
          }
     
          tower_by_tower_calib = recal_e[eta][phi];
      }
     
       const double raw_energy = raw_tower->get_energy();
       
       if(re_cal_flag == "CALIBMODE" && tower_by_tower_calib != 0)
       {
        recalibrated_e = raw_energy/tower_by_tower_calib;
       }
       else if(re_cal_flag == "CALIBMODE" && tower_by_tower_calib == 0)
       {  
         recalibrated_e = 0.;
       }
       else
       {
         recalibrated_e = raw_energy; 
       }
         
        const double calib_energy = (recalibrated_e - m_PedestalADC) * m_CalibConst_GeV_per_ADC;
 
        RawTower *calib_tower = new RawTowerv2(*raw_tower);
        calib_tower->set_energy(calib_energy);
        m_CalibTowerContainer->AddTower(key, calib_tower);
    }

    else
    {
      std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
                << " invalid calibration algorithm #" << m_CalibAlgorithmEnum
                << std::endl;
    gSystem->Exit(1);
    exit(1);
    }
  } 
 
  if (Verbosity())
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << "input sum energy = " << m_RawTowerContainer->getTotalEdep()
              << ", output sum digitalized value = "
              << m_CalibTowerContainer->getTotalEdep() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
  
}

void RawTowerCalibration::CreateNodes(PHCompositeNode *topNode)
{
  if (m_Detector.empty())
  {
    std::cout << PHWHERE << " Detector name not set, use RawTowerCalibration::Detector(std::string) to set" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << "Run Node missing, doing nothing." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  std::string TowerGeomNodeName = "TOWERGEOM_" + m_Detector;
  m_RawTowerGeomContainer = findNode::getClass<RawTowerGeomContainer>(topNode, TowerGeomNodeName);
  if (!m_RawTowerGeomContainer)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << " " << TowerGeomNodeName << " Node missing, doing bail out!"
              << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  if (Verbosity() >= 1)
  {
    m_RawTowerGeomContainer->identify();
  }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  std::string RawTowerNodeName = "TOWER_" + m_RawTowerNodePrefix + "_" + m_Detector;
  m_RawTowerContainer = findNode::getClass<RawTowerContainer>(dstNode, RawTowerNodeName);
  if (!m_RawTowerContainer)
  {
    std::cout << Name() << "::" << m_Detector << "::" << __PRETTY_FUNCTION__
              << " " << RawTowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", m_Detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(m_Detector);
    dstNode->addNode(DetNode);
  }

  // Be careful as a previous calibrator may have been registered for this detector
  std::string CalibTowerNodeName = "TOWER_" + m_CalibTowerNodePrefix + "_" + m_Detector;
  m_CalibTowerContainer = findNode::getClass<RawTowerContainer>(DetNode, CalibTowerNodeName);
  if (!m_CalibTowerContainer)
  {
    m_CalibTowerContainer = new RawTowerContainer(m_RawTowerContainer->getCalorimeterID());
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(m_CalibTowerContainer, CalibTowerNodeName, "PHObject");
    DetNode->addNode(towerNode);
  }
  return;
}
