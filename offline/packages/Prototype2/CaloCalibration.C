
#include "CaloCalibration.h"

#include "RawTower_Prototype2.h"
#include <g4cemc/RawTowerContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include "PROTOTYPE2_FEM.h"
#include <iostream>
#include <string>
#include <cassert>
#include "PROTOTYPE2_FEM.h"

using namespace std;

//____________________________________
CaloCalibration::CaloCalibration(const std::string& name) : //
    SubsysReco(string("CaloCalibration_") + name), //
    _calib_const_scale(-1), //negative pulse
    _calib_towers(NULL), _raw_towers(NULL), detector(name), //
    _calib_tower_node_prefix("CALIB"), //
    _raw_tower_node_prefix("RAW") //
{
}

//____________________________________
int
CaloCalibration::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int
CaloCalibration::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int
CaloCalibration::process_event(PHCompositeNode *topNode)
{

  if (verbosity)
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "Process event entered" << std::endl;
    }

  RawTowerContainer::ConstRange begin_end = _raw_towers->getTowers();
  RawTowerContainer::ConstIterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
    {
      const RawTowerDefs::keytype key = rtiter->first;
      const RawTower_Prototype2 *raw_tower =
          dynamic_cast<RawTower_Prototype2 *>(rtiter->second);
      assert(raw_tower);

      const double calibration_const = _calib_const_scale;

      vector<double> vec_signal_samples;
      for (int i = 0; i < RawTower_Prototype2::NSAMPLES; i++)
        {
          vec_signal_samples.push_back(
              calibration_const * raw_tower->get_signal_samples(i));
        }

      double peak = NAN;
      double peak_sample = NAN;
      double pedstal = NAN;

      PROTOTYPE2_FEM::SampleFit_PowerLawExp(vec_signal_samples, peak, peak_sample, pedstal);

      RawTower_Prototype2 *calib_tower = new RawTower_Prototype2(*raw_tower);
      calib_tower->set_energy(peak);
      calib_tower->set_time(peak_sample);

      for (int i = 0; i < RawTower_Prototype2::NSAMPLES; i++)
        {
          calib_tower->set_signal_samples(i, vec_signal_samples[i] - pedstal);
        }

      _calib_towers->AddTower(key, calib_tower);

    } //  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)

  if (verbosity)
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "input sum energy = " << _raw_towers->getTotalEdep()
          << ", output sum digitalized value = "
          << _calib_towers->getTotalEdep() << std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_______________________________________
void
CaloCalibration::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
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
      RawTowerNodeName.c_str());
  if (!_raw_towers)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " " << RawTowerNodeName << " Node missing, doing bail out!"
          << std::endl;
      throw std::runtime_error(
          "Failed to find " + RawTowerNodeName
              + " node in RawTowerCalibration::CreateNodes");
    }

  // Create the tower nodes on the tree
  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst(
      "PHCompositeNode", detector));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
    }

  // Be careful as a previous calibrator may have been registered for this detector
  CaliTowerNodeName = "TOWER_" + _calib_tower_node_prefix + "_" + detector;
  _calib_towers = findNode::getClass<RawTowerContainer>(DetNode,
      CaliTowerNodeName.c_str());
  if (!_calib_towers)
    {
      _calib_towers = new RawTowerContainer(
          RawTowerDefs::convert_name_to_caloid(detector));
      PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(
          _calib_towers, CaliTowerNodeName.c_str(), "PHObject");
      DetNode->addNode(towerNode);
    }
}

//___________________________________
int
CaloCalibration::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

