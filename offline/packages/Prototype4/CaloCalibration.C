#include "CaloCalibration.h"
#include "PROTOTYPE4_FEM.h"
#include "RawTower_Prototype4.h"

#include <TString.h>
#include <calobase/RawTowerContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

using namespace std;

//____________________________________
CaloCalibration::CaloCalibration(const std::string &name)
  : SubsysReco(string("CaloCalibration_") + name)
  , _calib_towers(NULL)
  , _raw_towers(NULL)
  , detector(name)
  , _calib_tower_node_prefix("CALIB")
  , _raw_tower_node_prefix("RAW")
  , _calib_params(name)
  , _fit_type(kPowerLawDoubleExpWithGlobalFitConstraint)
{
  SetDefaultParameters(_calib_params);
}

//____________________________________
int CaloCalibration::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________
int CaloCalibration::InitRun(PHCompositeNode *topNode)
{
  CreateNodeTree(topNode);

  if (Verbosity())
  {
    std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << " - print calibration parameters: " << endl;
    _calib_params.Print();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________
int CaloCalibration::process_event(PHCompositeNode *topNode)
{
  if (Verbosity())
  {
    std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << "Process event entered" << std::endl;
  }

  map<int, double> parameters_constraints;
  if (_fit_type == kPowerLawDoubleExpWithGlobalFitConstraint and _raw_towers->size() > 1)
  {
    if (Verbosity())
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
                << "Extract global fit parameter for constraining individual fits" << std::endl;
    }

    // signal template

    vector<double> vec_signal_samples(PROTOTYPE4_FEM::NSAMPLES, 0);

    int count = 0;

    RawTowerContainer::Range begin_end = _raw_towers->getTowers();
    RawTowerContainer::Iterator rtiter;
    for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
    {
      RawTower_Prototype4 *raw_tower =
          dynamic_cast<RawTower_Prototype4 *>(rtiter->second);
      assert(raw_tower);

      //      bool signal_check_pass = true;
      //      for (int i = 0; i < RawTower_Prototype4::NSAMPLES; i++)
      //      {
      //        if (raw_tower->get_signal_samples(i) <= 10 or raw_tower->get_signal_samples(i) >= ((1 << 14) - 10))
      //        {
      //          signal_check_pass = false;
      //          break;
      //        }
      //      }

      //      if (signal_check_pass)
      //      {
      ++count;

      for (int i = 0; i < RawTower_Prototype4::NSAMPLES; i++)
      {
        if (raw_tower->get_signal_samples(i) <= 10 or raw_tower->get_signal_samples(i) >= ((1 << 14) - 10))
          vec_signal_samples[i] = numeric_limits<double>::quiet_NaN();  // invalidate this sample
        else
          vec_signal_samples[i] += raw_tower->get_signal_samples(i);
      }
      //      }
    }

    if (count > 0)
    {
      for (int i = 0; i < RawTower_Prototype4::NSAMPLES; i++)
      {
        vec_signal_samples[i] /= count;
      }

      double peak = NAN;
      double peak_sample = NAN;
      double pedstal = NAN;
      map<int, double> parameters_io;

      PROTOTYPE4_FEM::SampleFit_PowerLawDoubleExp(vec_signal_samples, peak,
                                                  peak_sample, pedstal, parameters_io, Verbosity());
      //    std::map<int, double> &parameters_io,  //! IO for fullset of parameters. If a parameter exist and not an NAN, the fit parameter will be fixed to that value. The order of the parameters are
      //    ("Amplitude", "Sample Start", "Power", "Peak Time 1", "Pedestal", "Amplitude ratio", "Peak Time 2")

      parameters_constraints[1] = parameters_io[1];
      parameters_constraints[2] = parameters_io[2];
      parameters_constraints[3] = parameters_io[3];
      parameters_constraints[5] = parameters_io[5];
      parameters_constraints[6] = parameters_io[6];

      //      //special constraint if Peak Time 1 == Peak Time 2
      //      if (abs(parameters_constraints[6] - parameters_constraints[3]) < 0.1)
      //      {
      //        const double average_peak_time = (parameters_constraints[6] + parameters_constraints[3]) / 2.;
      //
      //        std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
      //                  << ": two shaping time are too close "
      //                  << parameters_constraints[3] << " VS " << parameters_constraints[6]
      //                  << ". Use average peak time instead: " << average_peak_time
      //                  << std::endl;
      //
      //        parameters_constraints[6] = average_peak_time;
      //        parameters_constraints[3] = average_peak_time;
      //        parameters_constraints[5] = 0;
      //      }
    }
    else
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
                << ": Failed to build signal template! Fit each channel individually instead" << std::endl;
    }
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
    RawTower_Prototype4 *raw_tower =
        dynamic_cast<RawTower_Prototype4 *>(rtiter->second);
    assert(raw_tower);

    double calibration_const = calib_const_scale;

    if (use_chan_calibration)
    {
      // channel to channel calibration.
      const int column = raw_tower->get_column();
      const int row = raw_tower->get_row();

      assert(column >= 0);
      assert(row >= 0);

      string calib_const_name(Form("calib_const_column%d_row%d", column, row));

      calibration_const *= _calib_params.get_double_param(calib_const_name);
    }

    vector<double> vec_signal_samples;
    for (int i = 0; i < RawTower_Prototype4::NSAMPLES; i++)
    {
      vec_signal_samples.push_back(
          raw_tower->get_signal_samples(i));
    }

    double peak = NAN;
    double peak_sample = NAN;
    double pedstal = NAN;

    switch (_fit_type)
    {
    case kPowerLawExp:
      PROTOTYPE4_FEM::SampleFit_PowerLawExp(vec_signal_samples, peak,
                                            peak_sample, pedstal, Verbosity());
      break;

    case kPeakSample:
      PROTOTYPE4_FEM::SampleFit_PeakSample(vec_signal_samples, peak,
                                           peak_sample, pedstal, Verbosity());
      break;

    case kPowerLawDoubleExp:
    {
      map<int, double> parameters_io;

      PROTOTYPE4_FEM::SampleFit_PowerLawDoubleExp(vec_signal_samples, peak,
                                                  peak_sample, pedstal, parameters_io, Verbosity());
    }
    break;

    case kPowerLawDoubleExpWithGlobalFitConstraint:
    {
      map<int, double> parameters_io(parameters_constraints);

      PROTOTYPE4_FEM::SampleFit_PowerLawDoubleExp(vec_signal_samples, peak,
                                                  peak_sample, pedstal, parameters_io, Verbosity());
    }
    break;
    default:
      cout << __PRETTY_FUNCTION__ << " - FATAL error - unkown fit type " << _fit_type << endl;
      exit(3);
      break;
    }

    // store the result - raw_tower
    if (std::isnan(raw_tower->get_energy()))
    {
      //Raw tower was never fit, store the current fit

      raw_tower->set_energy(peak);
      raw_tower->set_time(peak_sample);

      if (Verbosity())
      {
        raw_tower->identify();
      }
    }

    // store the result - calib_tower
    RawTower_Prototype4 *calib_tower = new RawTower_Prototype4(*raw_tower);
    calib_tower->set_energy(peak * calibration_const);
    calib_tower->set_time(peak_sample);

    for (int i = 0; i < RawTower_Prototype4::NSAMPLES; i++)
    {
      calib_tower->set_signal_samples(i, (vec_signal_samples[i] - pedstal) * calibration_const);
    }

    _calib_towers->AddTower(key, calib_tower);

  }  //  for (rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)

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
void CaloCalibration::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

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
                                                      RawTowerNodeName.c_str());
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
                                                        CaliTowerNodeName.c_str());
  if (!_calib_towers)
  {
    _calib_towers = new RawTowerContainer(_raw_towers->getCalorimeterID());
    PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(
        _calib_towers, CaliTowerNodeName.c_str(), "PHObject");
    DetNode->addNode(towerNode);
  }

  // update the parameters on the node tree
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  assert(parNode);
  const string paramnodename = string("Calibration_") + detector;

  //   this step is moved to after detector construction
  //   save updated persistant copy on node tree
  _calib_params.SaveToNodeTree(parNode, paramnodename);
}

//___________________________________
int CaloCalibration::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloCalibration::SetDefaultParameters(PHParameters &param)
{
  param.set_int_param("use_chan_calibration", 0);

  // additional scale for the calibration constant
  // negative pulse -> positive with -1
  param.set_double_param("calib_const_scale", 1);
}
