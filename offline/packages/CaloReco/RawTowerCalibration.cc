#include "RawTowerCalibration.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerv2.h>

#include <dbfile_calo_calib/CEmcCaloCalibSimpleCorrFilev1.h>
#include <dbfile_calo_calib/CaloCalibSimpleCorrFile.h>
#include <dbfile_calo_calib/HcalCaloCalibSimpleCorrFilev1.h>

#include <phparameter/PHParameters.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <cassert>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

RawTowerCalibration::RawTowerCalibration(const std::string &name)
  : SubsysReco(name)
  , _calib_algorithm(kNo_calibration)
  , detector("NONE")
  , _calib_tower_node_prefix("CALIB")
  , _raw_tower_node_prefix("RAW")
  , _tower_calib_params(name)
{
  //_tower_type = -1;
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
    std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
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

  if (_calib_algorithm == kDbfile_tbt_gain_corr)
  {
    if (detector.c_str()[0] == 'H')
      _cal_dbfile = (CaloCalibSimpleCorrFile *) new HcalCaloCalibSimpleCorrFilev1();
    else if (detector.c_str()[0] == 'C')
      _cal_dbfile = (CaloCalibSimpleCorrFile *) new CEmcCaloCalibSimpleCorrFilev1();
    else
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
                << "kDbfile_tbt_gain_corr  chosen but Detector Name not HCALOUT/IN or CEMC"
                << std::endl;
      return -999;
    }

    _cal_dbfile->Open(m_CalibrationFileName.c_str());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerCalibration::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity())
  {
    std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << "Process event entered" << std::endl;
  }

  RawTowerContainer::ConstRange begin_end = _raw_towers->getTowers();
  RawTowerContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
  {
    const RawTowerDefs::keytype key = rtiter->first;

    const RawTower *raw_tower = rtiter->second;
    assert(raw_tower);

    RawTowerGeom *raw_tower_geom = rawtowergeom->get_tower_geometry(
        raw_tower->get_id());
    assert(raw_tower_geom);

    if (_tower_type >= 0)
    {
      // Skip towers that don't match the type we are supposed to calibrate
      if (_tower_type != raw_tower_geom->get_tower_type())
      {
        continue;
      }
    }

    if (_calib_algorithm == kNo_calibration)
    {
      _calib_towers->AddTower(key, new RawTowerv2(*raw_tower));
    }
    else if (_calib_algorithm == kSimple_linear_calibration)
    {
      const double raw_energy = raw_tower->get_energy();
      const double calib_energy = (raw_energy - _pedstal_ADC) * _calib_const_GeV_ADC;

      RawTower *calib_tower = new RawTowerv2(*raw_tower);
      calib_tower->set_energy(calib_energy);
      _calib_towers->AddTower(key, calib_tower);
    }
    else if (_calib_algorithm == kTower_by_tower_calibration)
    {
      RawTowerDefs::CalorimeterId caloid = RawTowerDefs::decode_caloid(key);
      const int eta = raw_tower->get_bineta();
      const int phi = raw_tower->get_binphi();

      double tower_by_tower_calib = 1.;
      if (caloid == RawTowerDefs::LFHCAL)
      {
        const int l = raw_tower->get_binl();
        const std::string calib_const_name("calib_const_eta" + std::to_string(eta) + "_phi" + std::to_string(phi) + "_l" + std::to_string(l));

        tower_by_tower_calib = _tower_calib_params.get_double_param(calib_const_name);

        if (_pedestal_file == true)
        {
          const std::string pedstal_name("PedCentral_ADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi) + "_l" + std::to_string(l));
          _pedstal_ADC =
              _tower_calib_params.get_double_param(pedstal_name);
        }

        if (_GeV_ADC_file == true)
        {
          const std::string GeVperADCname("GeVperADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi) + "_l" + std::to_string(l));
          _calib_const_GeV_ADC =
              _tower_calib_params.get_double_param(GeVperADCname);
        }
      }
      else
      {
        const std::string calib_const_name("calib_const_eta" + std::to_string(eta) + "_phi" + std::to_string(phi));

        tower_by_tower_calib = _tower_calib_params.get_double_param(calib_const_name);

        if (_pedestal_file == true)
        {
          const std::string pedstal_name("PedCentral_ADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi));
          _pedstal_ADC =
              _tower_calib_params.get_double_param(pedstal_name);
        }

        if (_GeV_ADC_file == true)
        {
          const std::string GeVperADCname("GeVperADC_eta" + std::to_string(eta) + "_phi" + std::to_string(phi));
          _calib_const_GeV_ADC =
              _tower_calib_params.get_double_param(GeVperADCname);
        }
      }
      const double raw_energy = raw_tower->get_energy();
      const double calib_energy = (raw_energy - _pedstal_ADC) * _calib_const_GeV_ADC * tower_by_tower_calib;

      RawTower *calib_tower = new RawTowerv2(*raw_tower);
      calib_tower->set_energy(calib_energy);
      _calib_towers->AddTower(key, calib_tower);
    }
    //else if  // eventally this will be done exclusively of tow_by_tow
    else if (_calib_algorithm == kDbfile_tbt_gain_corr)
    {
      if (!_cal_dbfile)
      {
        std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
                  << "kDbfile_tbt_gain_corr  chosen but no file loaded" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      }

      float gain_factor = -888;
      //      gain_factor = _cal_dbfile->getCorr(key);

      const int eta = raw_tower->get_bineta();
      const int phi = raw_tower->get_binphi();

      gain_factor = _cal_dbfile->getCorr(eta, phi);

      const double raw_energy = raw_tower->get_energy();
      RawTower *calib_tower = new RawTowerv2(*raw_tower);

      // still include separate _calib_const_GeV_ADC factor
      // for global shifts.

      float corr_energy = raw_energy * gain_factor * _calib_const_GeV_ADC;
      calib_tower->set_energy(corr_energy);
      _calib_towers->AddTower(key, calib_tower);
    }
    else
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
                << " invalid calibration algorithm #" << _calib_algorithm
                << std::endl;

      return Fun4AllReturnCodes::ABORTRUN;
    }
  }  //  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)

  /*
  int towcount =0;
  RawTowerContainer::ConstRange begin_end2 = _calib_towers->getTowers();
  RawTowerContainer::ConstIterator rtiter2;
  std::cout << "Etowers  ---------------------================---------" 
	    << std::endl;

  for (rtiter2 = begin_end2.first; rtiter2 != begin_end2.second; ++rtiter2)
  {
    const RawTowerDefs::keytype key = rtiter2->first;
    const RawTower *aftcal_tower = rtiter2->second;

    if (towcount++ < 10)
      {
	std::cout << "E tow: " << key << "   " << aftcal_tower->get_energy()
                << std::endl;
      }


  }

  */

  if (Verbosity())
  {
    std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << "input sum energy = " << _raw_towers->getTotalEdep()
              << ", output sum digitalized value = "
              << _calib_towers->getTotalEdep() << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerCalibration::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawTowerCalibration::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << "Run Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find Run node in RawTowerCalibration::CreateNodes");
  }

  TowerGeomNodeName = "TOWERGEOM_" + detector;
  rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode,
                                                           TowerGeomNodeName);
  if (!rawtowergeom)
  {
    std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << " " << TowerGeomNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + TowerGeomNodeName + " node in RawTowerCalibration::CreateNodes");
  }

  if (Verbosity() >= 1)
  {
    rawtowergeom->identify();
  }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in RawTowerCalibration::CreateNodes");
  }

  RawTowerNodeName = "TOWER_" + _raw_tower_node_prefix + "_" + detector;
  _raw_towers = findNode::getClass<RawTowerContainer>(dstNode,
                                                      RawTowerNodeName);
  if (!_raw_towers)
  {
    std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << " " << RawTowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + RawTowerNodeName + " node in RawTowerCalibration::CreateNodes");
  }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst(
      "PHCompositeNode", detector));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(detector);
    dstNode->addNode(DetNode);
  }

  // Be careful as a previous calibrator may have been registered for this detector
  CaliTowerNodeName = "TOWER_" + _calib_tower_node_prefix + "_" + detector;
  _calib_towers = findNode::getClass<RawTowerContainer>(DetNode,
                                                        CaliTowerNodeName);
  if (!_calib_towers)
  {
    _calib_towers = new RawTowerContainer(_raw_towers->getCalorimeterID());
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_calib_towers, CaliTowerNodeName, "PHObject");
    DetNode->addNode(towerNode);
  }
  return;
}
