#include "CaloCalibration.h"
#include "PROTOTYPE2_FEM.h"
#include "RawTower_Prototype2.h"

#include <calobase/RawTowerContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/phool.h>
#include <phool/getClass.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <iostream>
#include <TString.h>
#include <cmath>
#include <string>
#include <cassert>
#include <cfloat>

using namespace std;

//____________________________________
CaloCalibration::CaloCalibration(const std::string& name) : //
    SubsysReco(string("CaloCalibration_") + name), //
    _calib_towers(NULL), _raw_towers(NULL), detector(name), //
    _calib_tower_node_prefix("CALIB"), //
    _raw_tower_node_prefix("RAW"), //
    _calib_params(name)
{
  SetDefaultParameters(_calib_params);
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

  if (Verbosity())
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " - print calibration parameters: "<<endl;
      _calib_params.Print();
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int
CaloCalibration::process_event(PHCompositeNode *topNode)
{

  if (Verbosity())
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "Process event entered" << std::endl;
    }

  const double calib_const_scale = _calib_params.get_double_param(
      "calib_const_scale");
  const bool use_chan_calibration = _calib_params.get_int_param(
      "use_chan_calibration") > 0;

  RawTowerContainer::Range begin_end = _raw_towers->getTowers();
  RawTowerContainer::Iterator rtiter;
  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
    {
      RawTowerDefs::keytype key = rtiter->first;
      RawTower_Prototype2 *raw_tower =
          dynamic_cast<RawTower_Prototype2 *>(rtiter->second);
      assert(raw_tower);

      double calibration_const = calib_const_scale;

      if (use_chan_calibration)
        {
          // channel to channel calibration.
          const int column = raw_tower->get_column();
          const int row = raw_tower->get_row();

          assert(column >= 0);
          assert(row >= 0);

          string calib_const_name(Form("calib_const_column%d_row%d",column,row));

          calibration_const *= _calib_params.get_double_param(calib_const_name);
        }

      vector<double> vec_signal_samples;
      for (int i = 0; i < RawTower_Prototype2::NSAMPLES; i++)
        {
          vec_signal_samples.push_back(
              raw_tower->get_signal_samples(i));
        }

      double peak = NAN;
      double peak_sample = NAN;
      double pedstal = NAN;

      PROTOTYPE2_FEM::SampleFit_PowerLawExp(vec_signal_samples, peak,
          peak_sample, pedstal, Verbosity());

      // store the result - raw_tower
      if (std::isnan(raw_tower->get_energy()))
        {
          //Raw tower was never fit, store the current fit

          raw_tower->set_energy(peak);
          raw_tower->set_time(peak_sample);
        }

      // store the result - calib_tower
      RawTower_Prototype2 *calib_tower = new RawTower_Prototype2(*raw_tower);
      calib_tower->set_energy(peak * calibration_const);
      calib_tower->set_time(peak_sample);

      for (int i = 0; i < RawTower_Prototype2::NSAMPLES; i++)
        {
          calib_tower->set_signal_samples(i, (vec_signal_samples[i] - pedstal) * calibration_const);
        }

      _calib_towers->AddTower(key, calib_tower);

    } //  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)

  if (Verbosity())
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
      _calib_towers = new RawTowerContainer(_raw_towers->getCalorimeterID());
      PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(
          _calib_towers, CaliTowerNodeName.c_str(), "PHObject");
      DetNode->addNode(towerNode);
    }

  // update the parameters on the node tree
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  assert(parNode);
  const string paramnodename = string("Calibration_") + detector;

  //   this step is moved to after detector construction
  //   save updated persistant copy on node tree
  _calib_params.SaveToNodeTree(parNode, paramnodename);

}

//___________________________________
int
CaloCalibration::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
CaloCalibration::SetDefaultParameters(PHParameters & param)
{

  param.set_int_param("use_chan_calibration", 0);

  // additional scale for the calibration constant
  // negative pulse -> positive with -1
  param.set_double_param("calib_const_scale", -1);

}
